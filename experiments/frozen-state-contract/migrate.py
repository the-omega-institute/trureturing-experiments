"""Project an old dual-write tree into the frozen state contract without repinning."""
import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import re
import sys
import tempfile

class MigrationError(ValueError):
    pass

HASH = re.compile(r'sha256:[0-9a-f]{64}\Z')
MODULE = re.compile(r'(?:[A-Z][A-Za-z0-9_]*[/])*[A-Z][A-Za-z0-9_]*\.lean\Z')

def ordinal(text):
    return text.encode('utf-16-be')

def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise MigrationError('duplicate JSON key: ' + key)
        result[key] = value
    return result

def read_json(path):
    try:
        return json.loads(path.read_bytes(), object_pairs_hook=unique_object)
    except (OSError, UnicodeError, ValueError) as error:
        raise MigrationError(str(path) + ': ' + str(error)) from error

def canonical(record):
    # System.Text.Json's default encoder: Basic Latin with HTML-sensitive characters escaped.
    def string(text):
        result = ['"']
        short = {'\b': r'\b', '\t': r'\t', '\n': r'\n', '\f': r'\f', '\r': r'\r', '\\': r'\\'}
        for char in text:
            if char in short:
                result.append(short[char])
            elif 0x20 <= ord(char) < 0x7f and char not in "\"&'<>+`":
                result.append(char)
            else:
                units = char.encode('utf-16-be')
                result.extend('\\u%04X' % int.from_bytes(units[i:i+2], 'big') for i in range(0, len(units), 2))
        result.append('"')
        return ''.join(result)
    def encode(value):
        if isinstance(value, str):
            return string(value)
        if isinstance(value, list):
            return '[' + ','.join(map(encode, value)) + ']'
        if isinstance(value, dict):
            return '{' + ','.join(string(key) + ':' + encode(item) for key, item in value.items()) + '}'
        raise MigrationError('unsupported canonical value')
    return (encode(record) + '\n').encode('ascii')

def project(payload):
    declarations = []
    seen = set()
    for declaration in payload['declaration_statement_ids']:
        if set(declaration) != {'declaration_name_key', 'kind', 'statement_id'}:
            raise MigrationError('invalid declaration fields')
        key, kind = declaration['declaration_name_key'], declaration['kind']
        if not isinstance(key, str) or not key or not isinstance(kind, str) or not kind:
            raise MigrationError('invalid declaration key or kind')
        if (key, kind) in seen:
            raise MigrationError('duplicate declaration: ' + key)
        seen.add((key, kind))
        declarations.append({'declaration_name_key': key, 'kind': kind})
    declarations.sort(key=lambda row: (ordinal(row['declaration_name_key']), ordinal(row['kind'])))
    return {'statement_id': payload['statement_id'], 'declarations': declarations}

def plan(root):
    root = Path(root).resolve()
    state_root = root / 'Golden/Frozen/state'
    accepted_root = root / 'Golden/Frozen/accepted'
    if not state_root.is_dir() or not accepted_root.is_dir():
        raise MigrationError('expected both old state and accepted directories')
    by_module = defaultdict(list)
    accepted_files = sorted(accepted_root.rglob('*'))
    for path in accepted_files:
        if path.is_dir():
            continue
        if path.parent != accepted_root or path.suffix != '.json' or path.is_symlink():
            raise MigrationError('unexpected accepted input: ' + str(path))
        event = read_json(path)
        if set(event) != {'event_hash', 'event_type', 'payload', 'schema_version'}:
            raise MigrationError('invalid event envelope: ' + str(path))
        if event['schema_version'] != 5 or event['event_type'] != 'Freeze':
            raise MigrationError('input is not a v5 Freeze: ' + str(path))
        payload = event['payload']
        if set(payload) != {'descriptor_selector', 'statement_id', 'declaration_statement_ids',
                            'prerequisite_frozen_node_ids'}:
            raise MigrationError('invalid Freeze payload: ' + str(path))
        by_module[payload['descriptor_selector']].append(payload)
    records, pins, count = {}, {}, 0
    for path in sorted(state_root.rglob('*')):
        if path.is_dir():
            continue
        selector = path.relative_to(state_root).as_posix()
        if not selector.endswith('.json') or path.is_symlink():
            raise MigrationError('unexpected state input: ' + str(path))
        selector = selector[:-5]
        if not MODULE.fullmatch(selector):
            raise MigrationError('invalid state module path: ' + selector)
        old = read_json(path)
        if set(old) != {'statement_id'} or not isinstance(old['statement_id'], str) or not HASH.fullmatch(old['statement_id']):
            raise MigrationError('expected old single-field state: ' + selector)
        matches = by_module[selector]
        if len(matches) != 1:
            raise MigrationError(f'{selector}: expected one Freeze, found {len(matches)}')
        if old['statement_id'] != matches[0]['statement_id']:
            raise MigrationError('state/Freeze pin mismatch: ' + selector)
        value = project(matches[0])
        records[path] = canonical(value)
        pins[selector] = old['statement_id']
        count += len(value['declarations'])
    # Accepted-only modules cannot be silently discarded by a whole-directory migration.
    if set(by_module) != set(pins):
        raise MigrationError('accepted/state member sets differ')
    migrated = {path.relative_to(state_root).as_posix()[:-5]: json.loads(wire)['statement_id']
                for path, wire in records.items()}
    if migrated != pins:
        raise MigrationError('projection changed members or pins')
    return records, {'members_before': len(pins), 'members_after': len(migrated),
                     'pins_equal': sum(pins[key] == migrated[key] for key in pins),
                     'declaration_entries': count, 'inconsistencies': 0}

def apply(root, records):
    # Recheck the entire source before replacing any member. Stage all bytes first.
    current, _ = plan(root)
    if current != records:
        raise MigrationError('inputs changed before application')
    backups = {path: path.read_bytes() for path in records}
    with tempfile.TemporaryDirectory(prefix='frozen-state-migration-') as temporary:
        stage = Path(temporary)
        for index, wire in enumerate(records.values()):
            (stage / str(index)).write_bytes(wire)
        applied = []
        try:
            for path, wire in records.items():
                with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
                    transient = Path(stream.name)
                    stream.write(wire)
                    stream.flush()
                    os.fsync(stream.fileno())
                os.replace(transient, path)
                applied.append(path)
        except BaseException:
            for path in reversed(applied):
                path.write_bytes(backups[path])
            raise
    # Directory retirement is explicit; no statement is computed from a current report.
    import shutil
    shutil.rmtree(Path(root) / 'Golden/Frozen/accepted')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repository', type=Path)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    try:
        records, result = plan(args.repository)
        if args.apply:
            apply(args.repository, records)
        print(json.dumps(result, sort_keys=True))
        return 0
    except MigrationError as error:
        print('FROZEN_STATE_MIGRATION_REJECTED ' + str(error), file=sys.stderr)
        return 1

if __name__ == '__main__':
    sys.exit(main())
