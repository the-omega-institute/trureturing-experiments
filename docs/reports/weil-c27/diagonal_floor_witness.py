"""A short actual c27 witness against the specified diagonal-floor debit.

Reuse the published directed scalar table, Arithmetic.even_entry and FLINT.
Only three missing middle scalars are generated. The negative number is an
upper bound for a sufficient Schur target, not a negative Weil test.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import zlib

import flint
from flint import acb, arb, ctx, fmpq

from arb_arithmetic import Arithmetic, grid_interval

HERE = Path(__file__).resolve().parent
SCALAR_SHA256 = 'f993be1aa48824cbf22276ac2d80e834ca9e54e4519e2040f437ebb4f724990f'
COEFFICIENTS = {1277: -1, 1278: 2, 1280: -1}
MIDDLE_MODES = range(1281, 1285)


def interval_ball(pair, bits):
    return arb(arb(fmpq(int(pair[0])+int(pair[1]), 1 << (bits+1))),
               arb(fmpq(int(pair[1])-int(pair[0]), 1 << (bits+1))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scalar-input', type=Path, default=HERE/'arithmetic-scalars.json.zlib')
    parser.add_argument('--precision', type=int, default=768)
    parser.add_argument('--output-bits', type=int, default=512)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output path.')
    raw = args.scalar_input.read_bytes()
    if hashlib.sha256(raw).hexdigest() != SCALAR_SHA256:
        raise ValueError('Use the published actual c27 scalar packet.')
    packet = json.loads(zlib.decompress(raw))
    bits = packet['grid_bits']
    if args.output_bits < 64 or args.precision < max(bits, args.output_bits)+64:
        parser.error('Use output-bits>=64 and at least 64 working guard bits.')
    supplier_hash = hashlib.sha256((HERE/'arb_arithmetic.py').read_bytes()).hexdigest()
    if supplier_hash != packet['source_sha256'] or packet['c'] != 27:
        raise ValueError('The existing scalar supplier does not match its published input.')
    ctx.prec = args.precision
    original = {row['n']: row for row in packet['rows']}
    arithmetic = Arithmetic(27, bits)
    used = sorted(set(COEFFICIENTS) | set(MIDDLE_MODES))
    data, new_rows, reused = {}, [], []
    for n in used:
        if n in original:
            row = original[n]
            reused.append(n)
        else:
            s, a = arithmetic.scalars(n)
            row = {'n': n, 's': grid_interval(s, bits), 'a': grid_interval(a, bits)}
            new_rows.append(row)
        data[n] = (interval_ball(row['s'], bits), interval_ball(row['a'], bits))

    retained = sum((arb(vi*vj)*arithmetic.even_entry(i, j, data)
                    for i, vi in COEFFICIENTS.items()
                    for j, vj in COEFFICIENTS.items()), arb(0))
    row_values = []
    for m in MIDDLE_MODES:
        coupling = sum((vi*arithmetic.even_entry(m, n, data)
                        for n, vi in COEFFICIENTS.items()), arb(0))
        # This is the existing G7 diagonal floor, sharper than the elementary
        # logarithmic lower comparison used by the archived far majorant.
        z = acb(arb(fmpq(1, 4)), arithmetic.pi*m/arithmetic.d)
        floor = z.digamma().real-arithmetic.pi.log()-arb(fmpq(59, 10))
        logarithmic = (arb(m)/arithmetic.d).log()-arithmetic.d/(arithmetic.pi*m)-arb(fmpq(59, 10))
        if not arb(fmpq(1, 16)) < floor < arb(fmpq(1, 15)):
            raise ArithmeticError('The four actual floor intervals do not fit the stated rational bounds.')
        if not 0 < logarithmic < arb(fmpq(1, 15)):
            raise ArithmeticError('The logarithmic comparison does not fit the rational upper bound.')
        row_values.append((m, coupling, floor, logarithmic))
    weakened_debit = 15*sum((coupling*coupling for _, coupling, _, _ in row_values), arb(0))
    upper_target = retained-weakened_debit
    if not retained > 30 or not upper_target < arb(fmpq(-13, 5)):
        raise ArithmeticError('The proposed rational vector does not verify the claimed failure margin.')
    export = lambda value: grid_interval(value, args.output_bits)
    result = {
        'schema': 'ACTUAL_C27_DIAGONAL_FLOOR_WITNESS_V1',
        'c': 27, 'last_retained': 1280,
        'basis': packet['basis'],
        'coefficients': [{'n': n, 'integer': v} for n, v in COEFFICIENTS.items()],
        'norm_squared': sum(v*v for v in COEFFICIENTS.values()),
        'middle_modes': list(MIDDLE_MODES),
        'floor_definition': 'Re psi(1/4+i*pi*m/log27)-log(pi)-59/10',
        'floor_lower': '1/16', 'floor_upper': '1/15',
        'weakened_inverse_weight': 15,
        'retained_quadratic': export(retained),
        'weakened_debit': export(weakened_debit),
        'sufficient_target_upper': export(upper_target),
        'strict_upper': '-13/5',
        'row_values': [{'m': m, 'coupling': export(coupling), 'floor': export(floor),
                        'logarithmic_floor': export(logarithmic)}
                       for m, coupling, floor, logarithmic in row_values],
        'scalar_grid_bits': bits, 'new_scalar_rows': new_rows, 'reused_modes': reused,
        'grid_bits': args.output_bits, 'working_bits': args.precision,
        'scalar_packet_sha256': SCALAR_SHA256, 'supplier_sha256': supplier_hash,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'python': platform.python_version(), 'python_flint': flint.__version__,
        'FLINT': flint.__FLINT_VERSION__,
        'scope': 'Conditional on the existing actual Fourier/form and G7 exterior '
                 'comparison. This witnesses failure of the N1280 diagonal-floor '
                 'Schur sufficient condition, including its weaker logarithmic '
                 'floor. The actual retained quadratic is positive. No negative '
                 'Weil test, full-form counterexample, Lean or RH certificate.'}
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print('retained', retained.str(32))
    print('weakened debit', weakened_debit.str(32))
    print('sufficient target upper', upper_target.str(32))
    print('result_sha256', hashlib.sha256(args.out.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()
