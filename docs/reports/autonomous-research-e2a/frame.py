#!/usr/bin/env python3
"""Immutable-source inventories only. Never evaluates a mathematical target."""
import argparse
import gzip
import json
from pathlib import Path
import re
import subprocess
import sys

BASE = '442e8d0a8855794547208a39a9ae17cd9798e572'


def git(repo, *args, data=None):
    return subprocess.run(['git', '-C', str(repo), *args], input=data,
                          stdout=subprocess.PIPE, check=True).stdout


def tree(repo):
    result = {}
    for row in git(repo, 'ls-tree', '-rz', BASE).split(b'\0'):
        if not row:
            continue
        header, path = row.split(b'\t', 1)
        mode, kind, oid = header.decode().split()
        if kind == 'blob':
            result[path.decode()] = oid
    return result


def blobs(repo, entries, paths):
    paths = sorted(paths)
    raw = git(repo, 'cat-file', '--batch',
              data=''.join(entries[p] + '\n' for p in paths).encode())
    offset = 0
    result = {}
    for path in paths:
        end = raw.index(b'\n', offset)
        oid, kind, length = raw[offset:end].split()
        if oid.decode() != entries[path] or kind != b'blob':
            raise ValueError('unexpected Git blob identity')
        start = end + 1
        size = int(length)
        result[path] = raw[start:start + size].decode('utf-8')
        offset = start + size + 1
    if offset != len(raw):
        raise ValueError('unconsumed Git bytes')
    return result


def output(value, path=None):
    data = (json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + '\n').encode()
    if path is None:
        sys.stdout.buffer.write(data)
    else:
        if path.suffix == '.gz':
            data = gzip.compress(data, mtime=0)
        path.write_bytes(data)


def decode_name(key):
    # The canonical key's component length is UTF-8 byte length.
    data = key.encode('utf-8')
    def parse(pos):
        if data[pos:pos+2] == b'n0':
            return [], pos + 2
        if data[pos:pos+3] not in (b'ns(', b'nn('):
            raise ValueError('unsupported canonical name form')
        numeric = data[pos:pos+3] == b'nn('
        parent, pos = parse(pos + 3)
        if data[pos:pos+1] != b',':
            raise ValueError('invalid canonical name parent')
        pos += 1
        if numeric:
            end = data.index(b')', pos)
            part = data[pos:end].decode()
        else:
            colon = data.index(b':', pos)
            size = int(data[pos:colon])
            pos = colon + 1
            end = pos + size
            part = data[pos:end].decode()
        if data[end:end+1] != b')':
            raise ValueError('invalid canonical name component length')
        return parent + [part], end + 1
    parts, end = parse(0)
    if end != len(data):
        raise ValueError('unconsumed canonical name')
    return '.'.join(parts)


def acquire(repo):
    locations = ['.lake', '.lake/build', '.lake/packages',
                 'tools/lean-inspector/.lake',
                 'tools/lean-inspector-interface/.lake',
                 'tools/lean-inspector-reg/.lake', 'Reg/.lake']
    status = {path: (repo / path).is_dir() for path in locations}
    return {'format': 'e2a-compiled-availability-v1', 'baseline': BASE,
            'directories_present': status,
            'complete_compiled_d5_value_graph': None,
            'complete_compiled_d5_value_graph_status': 'not-acquired',
            'reason': 'No compiled dependency graph is supplied; missing build '
                      'directories cannot be built in the read-only checkout.'}


def boundary(repo, entries, acquisition):
    if acquisition['baseline'] != BASE:
        raise ValueError('acquisition baseline mismatch')
    paths = [p for p in entries if p.startswith('Golden/Frozen/state/D5/')
             and p.endswith('.lean.json')]
    sources = blobs(repo, entries, paths)
    modules, rows, missing = [], [], []
    theorem_entries = private_entries = 0
    for path, source in sources.items():
        module_path = path.removeprefix('Golden/Frozen/state/').removesuffix('.json')
        module = module_path.removesuffix('.lean').replace('/', '.')
        state = json.loads(source)
        modules.append({'module': module, 'source': module_path,
                        'state': path, 'state_blob': entries[path],
                        'statement_id': state['statement_id'],
                        'source_present': module_path in entries})
        for item in state['declarations']:
            if item['kind'] != 'theorem':
                continue
            theorem_entries += 1
            key = item['declaration_name_key']
            try:
                name = decode_name(key)
            except (ValueError, UnicodeDecodeError) as exc:
                name = None
                missing.append({'source': path, 'declaration_name_key': key,
                                'reason': str(exc)})
            is_private = name is not None and name.startswith('_private.')
            if is_private:
                private_entries += 1
                continue
            row = {'declaration': name, 'module': module,
                   'visibility': 'public-name' if name else 'unresolved',
                   'users': None, 'user_count': None}
            if name is None:
                row['declaration_name_key'] = key
            rows.append(row)
    rows.sort(key=lambda r: (r['declaration'] or r['declaration_name_key'], r['module']))
    return {'format': 'e2a-boundary-frame-v1', 'baseline': BASE,
            'availability': acquisition,
            'population': {'frozen_d5_modules': len(modules),
                           'frozen_theorem_entries': theorem_entries,
                           'private_name_entries': private_entries,
                           'public_name_entries': sum(r['visibility'] == 'public-name' for r in rows),
                           'public_source_visibility': 'not-needed-after-acquisition-stop',
                           'verified_high_use_theorems': None,
                           'removable_hypothesis_candidates': None},
            'source_inventory_order': 'full declaration name ascending; module tie-break',
            'candidate_order': 'verified user count descending, full declaration name '
                               'ascending, original telescope index ascending',
            'candidate_order_status': 'unavailable-no-verified-counts',
            'all_listed_users_status': 'unmeasured',
            'hypothesis_expansion': 'not-reached-missing-proof-use',
            'modules': modules, 'theorems': rows, 'name_decode_missing': missing,
            'skipped': [], 'first_eligible': None, 'event': None,
            'status': 'infeasible-missing-declaration-proof-dependencies',
            'stop': 'No complete compiled D5 proof-value graph; no import or text proxy, '
                    'no ranking unknown as zero, no hypothesis expansion or replacement.'}



def mask_lean(source):
    """Preserve offsets/newlines, mask nested comments, strings and char literals."""
    chars = list(source)
    i, size = 0, len(source)
    while i < size:
        start = i
        if source.startswith('--', i):
            end = source.find('\n', i)
            i = size if end < 0 else end
        elif source.startswith('/-', i):
            depth = 1
            i += 2
            while i < size and depth:
                if source.startswith('/-', i):
                    depth += 1
                    i += 2
                elif source.startswith('-/', i):
                    depth -= 1
                    i += 2
                else:
                    i += 1
            if depth:
                raise ValueError('unterminated Lean comment')
        elif source[i] == '"':
            i += 1
            while i < size:
                if source[i] == '\\':
                    i += 2
                elif source[i] == '"':
                    i += 1
                    break
                else:
                    i += 1
            else:
                raise ValueError('unterminated Lean string')
        elif (source[i] == "'" and
              (i == 0 or not (source[i-1].isalnum() or source[i-1] in "_'")) and
              re.match(r"'(?:[^'\\\n]|\\.)'", source[i:])):
            i += len(re.match(r"'(?:[^'\\\n]|\\.)'", source[i:])[0])
        else:
            i += 1
            continue
        for j in range(start, min(i, size)):
            if chars[j] != '\n':
                chars[j] = ' '
    return ''.join(chars)


IDENT = r'(?:«[^»\n]+»|[\w][\w\'!?₀-₉]*)(?:\.(?:«[^»\n]+»|[\w][\w\'!?₀-₉]*))*'
COMMAND = re.compile(
    r'(?m)^[ \t]*(?:(?:private|protected|noncomputable|unsafe|partial|public|local)\s+)*'
    r'(namespace|section|end|def|abbrev|theorem|lemma|instance|opaque|axiom|structure|class|inductive)'
    r'\b(?:[ \t]+(' + IDENT + r'))?')
UNKNOWN = re.compile(r'\bcontinuation\s*:=\s*\.\s*unknown\b')
ANY_CONTINUATION = re.compile(r'\bcontinuation\s*:=')


def residual(repo, entries):
    paths = sorted(p for p in entries if p.endswith('.lean')
                   and (p.startswith('Reg/') or p.startswith('D5/')))
    rows, missing, fields, files_with_fields = [], [], 0, 0
    all_field_locations = []
    for offset in range(0, len(paths), 200):
        sources = blobs(repo, entries, paths[offset:offset+200])
        for path, source in sources.items():
            if 'continuation' not in source:
                continue
            masked = mask_lean(source)
            all_fields = list(ANY_CONTINUATION.finditer(masked))
            fields += len(all_fields)
            files_with_fields += bool(all_fields)
            all_field_locations.extend({'source': path, 'line': masked.count('\n', 0, x.start())+1,
                                        'explicit_unknown': bool(UNKNOWN.match(masked, x.start()))}
                                       for x in all_fields)
            occurrences = list(UNKNOWN.finditer(masked))
            if not occurrences:
                continue
            commands = list(COMMAND.finditer(masked))
            scope = []
            declarations = []
            for index, command in enumerate(commands):
                kind, name = command.group(1), command.group(2)
                if kind == 'namespace':
                    if not name:
                        missing.append({'source': path, 'line': masked.count('\n', 0, command.start())+1,
                                        'reason': 'namespace-name-unresolved'})
                    scope.append(('namespace', name or 'UNRESOLVED'))
                elif kind == 'section':
                    scope.append(('section', name))
                elif kind == 'end':
                    if scope:
                        scope.pop()
                    else:
                        missing.append({'source': path, 'line': masked.count('\n', 0, command.start())+1,
                                        'reason': 'unmatched-scope-end'})
                else:
                    namespaces = [n for k, n in scope if k == 'namespace']
                    if name and name.startswith('_root_.'):
                        full_name = name[len('_root_.'):]
                    else:
                        full_name = '.'.join(namespaces + [name or 'UNRESOLVED'])
                    end = commands[index+1].start() if index+1 < len(commands) else len(source)
                    declarations.append((command.start(), end, full_name, kind, name))
            for occurrence in occurrences:
                owners = [d for d in declarations if d[0] <= occurrence.start() < d[1]]
                line = masked.count('\n', 0, occurrence.start())+1
                if len(owners) != 1:
                    missing.append({'source': path, 'line': line,
                                    'reason': 'unknown-field-owner-unresolved'})
                    continue
                start, end, name, kind, local = owners[0]
                declaration_text = source[start:end].strip()
                declaration_mask = masked[start:end]
                header_end = re.search(r':=|\bwhere\b', declaration_mask)
                header = declaration_text[:header_end.start()].strip() if header_end else declaration_text
                is_registration = bool(re.search(r'\bRegistration\b', header))
                if kind not in ('def', 'abbrev') or not local or not is_registration:
                    missing.append({'source': path, 'line': line, 'declaration': name,
                                    'reason': 'field-owner-not-explicitly-typed-registration'})
                row = {'registration': name, 'declaration_kind': kind,
                       'source': path, 'source_blob': entries[path],
                       'declaration_line': masked.count('\n', 0, start)+1,
                       'continuation_line': line,
                       'continuation': '.unknown', 'declaration_header': header,
                       'registration_type_syntactic': is_registration,
                       'assessment': 'not-assessed-after-prefix', 'skip': None}
                # Preserve field text without extracting/probing proof bodies.
                for field in ['unitName', 'realizationName', 'realizationSource', 'arena',
                              'objectArena', 'catalog', 'localNames', 'realization',
                              'readout', 'escapeFrom', 'familyRecord']:
                    match = re.search(r'(?m)^[ \t]*' + field + r'\s*:=([^\n]*)', declaration_mask)
                    row[field] = source[start+match.start(1):start+match.end(1)].strip() if match else None
                rows.append(row)
    rows.sort(key=lambda r: (r['registration'], r['source'], r['continuation_line']))
    for index, row in enumerate(rows, 1):
        row['rank'] = index
    duplicates = {}
    for row in rows:
        duplicates.setdefault(row['registration'], []).append(row)
    collisions = [{'registration': name,
                   'occurrences': [{'source': r['source'], 'line': r['continuation_line']}
                                   for r in group]}
                  for name, group in sorted(duplicates.items()) if len(group) > 1]
    return {'format': 'e2a-residual-inventory-v1', 'baseline': BASE,
            'population': {'lean_files_scanned': len(paths),
                           'files_with_explicit_continuation_fields': files_with_fields,
                           'explicit_continuation_fields': fields,
                           'explicit_unknown_registration_rows': len(rows),
                           'owner_parse_missing': len(missing)},
            'scope': 'all Git-tracked Reg/**/*.lean and D5/**/*.lean at baseline',
            'mechanism': 'offset-preserving comment/string mask, explicit field tokens, '
                         'source namespace/section/declaration ownership; not compiled validation',
            'order': 'full registration declaration name ascending; path and line expose collisions',
            'all_explicit_continuation_fields': all_field_locations,
            'registrations': rows, 'parse_missing': missing, 'name_collisions': collisions,
            'first_eligible': None, 'event': None,
            'status': 'inventory-only-awaiting-bounded-event-assessment'}


def apply_decisions(inventory, decisions):
    if decisions['baseline'] != BASE:
        raise ValueError('decision baseline mismatch')
    assessments = decisions['assessments']
    expected_ranks = list(range(1, len(assessments) + 1))
    if [a['rank'] for a in assessments] != expected_ranks or len(assessments) > 20:
        raise ValueError('non-prefix or over-budget assessment')
    if inventory['parse_missing'] or inventory['name_collisions']:
        raise ValueError('source ownership incomplete')
    selected = None
    for item in assessments:
        row = inventory['registrations'][item['rank'] - 1]
        for key in ['registration', 'source', 'source_blob']:
            if row[key] != item[key]:
                raise ValueError('assessment source identity mismatch: ' + key)
        row['assessment'] = item['status']
        row['assessment_evidence'] = item
        if selected is not None:
            raise ValueError('assessment after first eligible candidate')
        if item['status'] == 'eligible-frozen':
            selected = row['registration']
        elif item['status'] == 'skipped-known-result':
            row['skip'] = item['reason']
        else:
            raise ValueError('unsupported assessment status')
    unit = decisions['slots']['residual']
    if selected != unit['source_registration']:
        raise ValueError('selected event source differs from first eligible row')
    inventory.update(format='e2a-residual-frame-v1', first_eligible=selected,
                     event=unit['event'],
                     status=unit['status'], assessed_prefix=len(assessments),
                     skipped_count=sum(a['status'].startswith('skipped-') for a in assessments),
                     remaining_status='not-event-assessed-after-first-eligible')
    return inventory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--main-repo', required=True, type=Path)
    parser.add_argument('mode', choices=['acquire', 'boundary', 'residual'])
    parser.add_argument('--availability', type=Path)
    parser.add_argument('--decisions', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    entries = tree(args.main_repo)
    if args.mode == 'acquire':
        output(acquire(args.main_repo), args.output)
    elif args.mode == 'residual':
        inventory = residual(args.main_repo, entries)
        if args.decisions:
            inventory = apply_decisions(inventory, json.loads(args.decisions.read_text()))
        output(inventory, args.output)
    elif args.mode == 'boundary':
        if args.availability is None:
            parser.error('--availability is required for boundary')
        output(boundary(args.main_repo, entries,
                        json.loads(args.availability.read_text())), args.output)


if __name__ == '__main__':
    main()
