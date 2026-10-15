"""Check c9 source moments and replay the hash-bound retained comparison."""
import argparse
import base64
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
import zlib
from pathlib import Path


def check_bytes(payload, record, label):
    if len(payload) != record['bytes']:
        raise ValueError(f'Wrong byte count: {label}')
    if hashlib.sha256(payload).hexdigest() != record['sha256']:
        raise ValueError(f'Wrong SHA-256: {label}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--moment-checker', type=Path,
                        help='Pinned external source-only moment oracle; required for full replay.')
    parser.add_argument('--check-inputs-only', action='store_true',
                        help='Check package bindings without a sign calculation.')
    args = parser.parse_args()
    package = Path(__file__).resolve().parent
    manifest_path = package / 'input-manifest.json'
    manifest = json.loads(manifest_path.read_bytes())
    if manifest['schema'] != 'ACTUAL_C9_RETAINED_INPUTS_V2':
        raise ValueError('Wrong scientific manifest.')
    versions = {'python': platform.python_version(),
                'python_flint': importlib.metadata.version('python-flint'),
                'numpy': importlib.metadata.version('numpy')}
    if versions != manifest['arithmetic_versions']:
        raise ValueError(f'Expected {manifest["arithmetic_versions"]}; got {versions}')
    for name, record in manifest['programs'].items():
        check_bytes((package / name).read_bytes(), record, name)
    check_bytes((package / manifest['consumer']['path']).read_bytes(),
                manifest['consumer'], 'consumer')
    decoded = {}
    for key in ['prime', 'kernel', 'moments']:
        record = manifest[key]
        packed = record['packed']
        payload = (package / packed['path']).read_bytes()
        check_bytes(payload, packed, packed['path'])
        compressed = base64.b64decode(b''.join(payload.splitlines()), validate=True)
        decoded[key] = zlib.decompress(compressed)
        check_bytes(decoded[key], record, key)
    if args.check_inputs_only:
        print('Package bindings verified; no sign calculation executed.')
        return
    if args.moment_checker is None:
        parser.error('--moment-checker is required for full replay')
    destination = args.out_dir.resolve()
    if destination == package:
        raise ValueError('Use a result directory outside the published package.')
    destination.mkdir(parents=True, exist_ok=True)
    for key, payload in decoded.items():
        (destination / manifest[key]['path']).write_bytes(payload)
    fresh_moments = destination / 'fresh-source-moments.json'
    subprocess.run([sys.executable, str(package / 'c9_source_moments.py'),
                    '--checker', str(args.moment_checker), '--out', str(fresh_moments)],
                   check=True)
    fresh = json.loads(fresh_moments.read_bytes())
    retained = json.loads(decoded['moments'])
    for key in ['parameters', 'source_checker', 'rows']:
        if fresh[key] != retained[key]:
            raise ValueError(f'Regenerated moment input differs: {key}')
    command = [sys.executable, str(package / manifest['consumer']['path']),
               '--prime', str(destination / manifest['prime']['path']),
               '--kernel', str(destination / manifest['kernel']['path']),
               '--moments', str(destination / manifest['moments']['path']),
               '--manifest', str(manifest_path),
               '--out', str(destination / 'result.json'), '--bits', '1536']
    subprocess.run(command, check=True)
    print(f'Full retained trial and directed pivots: {destination / "result.json"}')


if __name__ == '__main__':
    main()
