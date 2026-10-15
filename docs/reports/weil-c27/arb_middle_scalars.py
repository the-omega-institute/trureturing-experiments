"""Supply only missing actual c27 middle scalars through the existing oracle.

The published low packet and three previously supplied middle rows are reused.
No arithmetic formula or special-function algorithm is reimplemented.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time
import zlib

import flint
from flint import ctx

from arb_arithmetic import Arithmetic, grid_interval
from arb_far_coupling import SCALAR_SHA256, pinned_packet


HERE = Path(__file__).resolve().parent
WITNESS_SHA256 = 'f45f3a333b50aa60c38c5b2ee9be416a165a1c719ba7dfb3ca6dfbe251bf1840'


def published_rows(scalar_path, witness_path):
    scalar = pinned_packet(scalar_path, SCALAR_SHA256)
    raw = witness_path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != WITNESS_SHA256:
        raise ValueError('Use the published four-row witness.')
    witness = json.loads(raw)
    supplier = hashlib.sha256((HERE/'arb_arithmetic.py').read_bytes()).hexdigest()
    if (scalar['c'] != 27 or witness['c'] != 27
            or scalar['source_sha256'] != supplier
            or witness['supplier_sha256'] != supplier
            or witness['scalar_packet_sha256'] != SCALAR_SHA256
            or scalar['grid_bits'] != witness['scalar_grid_bits']):
        raise ValueError('Published sources or scalar grids disagree.')
    rows = {row['n']: row for row in scalar['rows']}
    for row in witness['new_scalar_rows']:
        if row['n'] in rows and rows[row['n']] != row:
            raise ValueError('Conflicting published row.')
        rows[row['n']] = row
    return scalar, rows


def read_middle(path, sha256, supplier_hash, bits, precision):
    packet = pinned_packet(path, sha256)
    if (packet['schema'] != 'ACTUAL_C27_MIDDLE_SCALARS_V1'
            or packet['c'] != 27 or packet['first_middle'] != 1281
            or packet['grid_bits'] != bits
            or packet['working_bits'] != precision
            or packet['supplier_sha256'] != supplier_hash
            or packet['source_sha256'] != hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            or packet['scalar_packet_sha256'] != SCALAR_SHA256
            or packet['witness_sha256'] != WITNESS_SHA256):
        raise ValueError('Middle cache belongs to other parameters or sources.')
    rows = {}
    for row in packet['rows']:
        n = row['n']
        if type(n) is not int or not 1281 <= n <= packet['last_middle'] or n in rows:
            raise ValueError('Invalid or duplicate cached mode.')
        for family in ('s', 'a'):
            pair = row[family]
            if (len(pair) != 2 or any(type(v) is not str or str(int(v)) != v for v in pair)
                    or not 0 <= int(pair[1])-int(pair[0]) <= 2):
                raise ValueError('Invalid cached scalar interval.')
        rows[n] = row
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scalar-input', type=Path, default=HERE/'arithmetic-scalars.json.zlib')
    parser.add_argument('--witness-input', type=Path, default=HERE/'diagonal-floor-witness.json')
    parser.add_argument('--last-middle', type=int, default=32768)
    parser.add_argument('--precision', type=int, default=768)
    parser.add_argument('--reuse-input', type=Path)
    parser.add_argument('--reuse-sha256')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output path.')
    scalar, original = published_rows(args.scalar_input, args.witness_input)
    bits = scalar['grid_bits']
    if not 1281 <= args.last_middle <= 32768 or args.precision < bits+64:
        parser.error('Require1281<=last-middle<=32768 and at least64 guard bits.')
    if bool(args.reuse_input) != bool(args.reuse_sha256):
        parser.error('Pin --reuse-input with --reuse-sha256.')
    ctx.prec = args.precision
    supplier = hashlib.sha256((HERE/'arb_arithmetic.py').read_bytes()).hexdigest()
    reused = (read_middle(args.reuse_input, args.reuse_sha256, supplier, bits, args.precision)
              if args.reuse_input else {})
    if set(reused) & set(original):
        raise ValueError('The middle cache must contain only newly supplied rows.')
    if any(n > args.last_middle for n in reused):
        raise ValueError('The reused cache exceeds the requested range.')
    started = time.monotonic()
    arithmetic = Arithmetic(27, bits)
    rows, generated = [], 0
    for n in range(1281, args.last_middle+1):
        if n in original:
            continue
        if n in reused:
            row = reused[n]
        else:
            s, a = arithmetic.scalars(n)
            row = {'n': n, 's': grid_interval(s, bits), 'a': grid_interval(a, bits)}
            generated += 1
        rows.append(row)
        if generated and generated % 512 == 0 and n not in reused:
            print('generated', generated, 'through', n, 'seconds', time.monotonic()-started, flush=True)
    packet = {
        'schema': 'ACTUAL_C27_MIDDLE_SCALARS_V1', 'c': 27,
        'first_middle': 1281, 'last_middle': args.last_middle,
        'grid_bits': bits, 'working_bits': args.precision,
        'rows': rows, 'row_count': len(rows), 'generated_rows': generated,
        'reused_cache_rows': len(reused),
        'reused_published_modes': sorted(n for n in original if 1281 <= n <= args.last_middle),
        'scalar_packet_sha256': SCALAR_SHA256, 'witness_sha256': WITNESS_SHA256,
        'supplier_sha256': supplier,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'python': platform.python_version(), 'python_flint': flint.__version__,
        'FLINT': flint.__FLINT_VERSION__,
        'elapsed_seconds': time.monotonic()-started,
        'scope': 'New actual middle scalar inputs only. Existing Fourier/form '
                 'identification is an input; no matrix sign, full residual '
                 'consumer, Lean or RH certificate.'}
    expected = set(range(1281, args.last_middle+1))-set(original)
    if {row['n'] for row in rows} != expected:
        raise ArithmeticError('Incomplete actual middle scalar coverage.')
    raw = (json.dumps(packet, sort_keys=True, separators=(',', ':'))+'\n').encode()
    args.out.write_bytes(zlib.compress(raw, level=9))
    print('rows', len(rows), 'new', generated, 'seconds', packet['elapsed_seconds'])
    print('packet_sha256', hashlib.sha256(args.out.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()
