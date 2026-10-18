#!/usr/bin/env python3
"""Check the bytes, JSON structure and finite counts of the certificate data."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parent


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'duplicate JSON key: {key}')
        result[key] = value
    return result


def reject_constant(value):
    raise ValueError(f'invalid JSON constant: {value}')


def read_json(raw):
    return json.loads(raw, object_pairs_hook=unique_object, parse_constant=reject_constant)


def read_payload(raw, format_name):
    if format_name == 'jsonl':
        value = [read_json(line) for line in raw.splitlines() if line.strip()]
        return value, len(value)
    if format_name != 'json':
        raise ValueError(f'unknown payload format: {format_name}')
    value = read_json(raw)
    return value, len(value) if isinstance(value, list) else 1


def explicit_counts(value):
    result = Counter(root_boxes=0, arcs=0, leaves=0)

    def visit(node):
        if isinstance(node, dict):
            if isinstance(node.get('center'), list) and isinstance(node.get('R'), list):
                result['root_boxes'] += 1
            for key, child in node.items():
                if key in ('arcs', 'leaves') and isinstance(child, list):
                    result[key] += len(child)
                visit(child)
        elif isinstance(node, list):
            for child in node:
                visit(child)

    visit(value)
    return dict(result)


def check_dataset(directory, source_repo):
    manifest = read_json((directory / 'manifest.json').read_bytes())
    entries = manifest['entries']
    if len(entries) != manifest['payload_files']:
        raise ValueError(f'{directory.name}: payload file count mismatch')
    if sum(entry['bytes'] for entry in entries) != manifest['payload_bytes']:
        raise ValueError(f'{directory.name}: total payload byte count mismatch')
    expected = {entry['path'] for entry in entries}
    actual = {
        path.relative_to(directory).as_posix()
        for path in directory.rglob('*')
        if path.is_file() and path.name != 'manifest.json'
    }
    if len(expected) != len(entries) or actual != expected:
        raise ValueError(f'{directory.name}: missing, duplicate or unlisted payload')

    source = None
    if source_repo is not None:
        source = subprocess.check_output([
            'git', '-C', str(source_repo), 'show',
            manifest['source_commit'] + ':' + manifest['source_path'],
        ])
        if len(source) != manifest['source_bytes'] or hashlib.sha256(source).hexdigest() != manifest['source_sha256']:
            raise ValueError(f'{directory.name}: source byte identity mismatch')

    totals = Counter(files=0, bytes=0, records=0, root_boxes=0, arcs=0, leaves=0)
    for entry in entries:
        payload = (directory / entry['path']).read_bytes()
        if len(payload) != entry['bytes'] or hashlib.sha256(payload).hexdigest() != entry['sha256']:
            raise ValueError(f"{directory.name}/{entry['path']}: byte identity mismatch")
        if source is not None:
            start, end = entry['source_body_byte_range']
            if source[start:end] != payload:
                raise ValueError(f"{directory.name}/{entry['path']}: original fence mismatch")
        value, records = read_payload(payload, entry['format'])
        counts = explicit_counts(value)
        if records != entry['records'] or any(counts[key] != entry[key] for key in counts):
            raise ValueError(f"{directory.name}/{entry['path']}: finite record count mismatch")
        totals.update(files=1, bytes=len(payload), records=records, **counts)
    print(json.dumps({'dataset': directory.name, **totals}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-repo', type=Path)
    args = parser.parse_args()
    for name in ('connected-interval', 'two-parameter-strip'):
        check_dataset(ROOT / name, args.source_repo)


if __name__ == '__main__':
    main()
