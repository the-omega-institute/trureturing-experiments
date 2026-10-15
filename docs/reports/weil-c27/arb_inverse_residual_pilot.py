"""Test the existing inverse-residual bound on the known c27 debit witness.

All middle rows and the infinite far Gram belong to the same actual columns.
The old signed moments and arithmetic oracle are reused. This tests an
estimate, not previously unknown positivity of the selected high vector.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
from types import SimpleNamespace
import zlib

import flint
from flint import arb, arb_mat, ctx, fmpq

from arb_arithmetic import Arithmetic, grid_interval
from arb_far_coupling import far_bounds, pinned_packet, reconstruct_factors, upper_dyadic
from arb_middle_scalars import HERE, SCALAR_SHA256, WITNESS_SHA256, published_rows, read_middle
from diagonal_floor_witness import COEFFICIENTS, interval_ball


MIDDLE_SHA256 = '55212d0181673c015de503928265509e594eb0f81b85e3afb34cc28a103d18db'
FAR_SHA256 = 'ac35edacd7f55de470586e3b8a5915bcaa49b9d2c4c9c2ee073eeaaa8c9c73e7'
EXTERIOR = (1281, 1282, 1283, 1284)


def trial_integer(value, bits):
    mantissa, exponent = map(int, value.mid().man_exp())
    exponent += bits
    return (mantissa << exponent if exponent >= 0 else
            (mantissa+(1 << (-exponent-1))) >> (-exponent))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scalar-input', type=Path, default=HERE/'arithmetic-scalars.json.zlib')
    parser.add_argument('--witness-input', type=Path, default=HERE/'diagonal-floor-witness.json')
    parser.add_argument('--middle-input', type=Path, default=HERE/'middle-scalars.json.zlib')
    parser.add_argument('--far-input', type=Path, default=HERE/'far-coupling.json.zlib')
    parser.add_argument('--precision', type=int, default=768)
    parser.add_argument('--output-bits', type=int, default=512)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output path.')
    if args.output_bits < 64 or args.precision < max(640, args.output_bits)+64:
        parser.error('Require output-bits>=64 and at least64 working guard bits.')
    ctx.prec = args.precision
    started = time.monotonic()
    scalar, rows = published_rows(args.scalar_input, args.witness_input)
    middle = pinned_packet(args.middle_input, MIDDLE_SHA256)
    new_rows = read_middle(args.middle_input, MIDDLE_SHA256,
                          scalar['source_sha256'], scalar['grid_bits'], middle['working_bits'])
    if set(new_rows) & set(rows) or middle['last_middle'] != 32768:
        raise ValueError('Wrong middle coverage or duplicate published modes.')
    rows.update(new_rows)
    far = pinned_packet(args.far_input, FAR_SHA256)
    if (far['N'], far['M'], far['jet_order'], far['moment_diagonals'],
            far['grid_bits']) != (1280, 32768, 56, 56, 640):
        raise ValueError('Wrong published far-moment parameters.')
    for key, name in [('source_sha256', 'arb_far_coupling.py'),
                      ('moment_supplier_sha256', 'arb_tail_moments.py')]:
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != far[key]:
            raise ValueError('The far supplier identity changed.')
    if not set(range(32769)) <= set(rows):
        raise ValueError('Missing actual low or middle scalar input.')
    bits = scalar['grid_bits']
    data = {n: (interval_ball(row['s'], bits), interval_ball(row['a'], bits))
            for n, row in rows.items()}
    arithmetic = Arithmetic(27, bits)
    v = COEFFICIENTS
    selected = EXTERIOR
    size = len(selected)+1
    retained = sum((vi*vj*arithmetic.even_entry(i, j, data)
                    for i, vi in v.items() for j, vj in v.items()), arb(0))
    B = [sum((vi*arithmetic.even_entry(j, n, data) for n, vi in v.items()), arb(0))
         for j in selected]
    F = [[arithmetic.even_entry(i, j, data) for j in selected] for i in selected]

    # Reuse the same168 moment intervals. Extend only the feature columns and
    # recompute all geometric and arithmetic amplification allowances.
    extended = dict(far, N=max(selected))
    extended_scalar = dict(scalar, rows=[rows[n] for n in range(max(selected)+1)])
    alpha, D, P, H, _old_error = reconstruct_factors(extended, extended_scalar)
    features = [[sum((vi*P[k, n] for n, vi in v.items()), fmpq(0))]
                +[P[k, n] for n in selected] for k in range(P.nrows())]
    feature_matrix = flint.fmpq_mat(features)
    far_center = D*alpha*(feature_matrix.transpose()*H*feature_matrix)
    scalar_radius = arb(max(int(rows[n]['s'][1])-int(rows[n]['s'][0])
                            for n in range(max(selected)+1)))/(1 << (bits+1))
    bounds_input = SimpleNamespace(M=far['M'], L=arithmetic.d, pi=arithmetic.pi,
                                  budget=2*arithmetic.H+sum((w for w, _ in arithmetic.primes), arb(0)))
    bounds = far_bounds(bounds_input, max(selected), far['jet_order'], scalar_radius,
                        far['grid_bits'], far['moment_diagonals'])
    # D is the archived exact dyadic far weight, not a newly selected midpoint.
    far_error = arb(D)*bounds['unweighted_matrix_error']
    norms = [arb(sum(c*c for c in v.values())).sqrt()]+[arb(1)]*len(selected)
    gram = [[arb(far_center[i, j])+arb(0, (far_error*norms[i]*norms[j]).upper())
             for j in range(size)] for i in range(size)]
    weights = []
    weight_bits = 256
    for m in range(1281, far['M']+1):
        h = (arb(m)/arithmetic.d).log()-arithmetic.d/(arithmetic.pi*m)-arb(fmpq(59, 10))
        if not h > 0:
            raise ArithmeticError('Nonpositive actual middle logarithmic floor.')
        numerator = upper_dyadic(1/h, weight_bits)
        weight = arb(fmpq(int(numerator), 1 << weight_bits))
        weights.append(numerator)
        column = [sum((vi*arithmetic.even_entry(m, n, data)
                       for n, vi in v.items()), arb(0))]
        column += [arithmetic.even_entry(m, n, data) for n in selected]
        for i in range(size):
            for j in range(i+1):
                gram[i][j] += weight*column[i]*column[j]
        if m % 4096 == 0:
            print('middle through', m, 'seconds', time.monotonic()-started, flush=True)
    for i in range(size):
        for j in range(i):
            gram[j][i] = gram[i][j]

    cases = []
    for count in (1, 2, 4):
        K = [[gram[i+1][j+1]-F[i][j] for j in range(count)] for i in range(count)]
        T = [gram[i+1][0]-B[i] for i in range(count)]
        proposal = arb_mat(K).solve(arb_mat([[x] for x in T]))
        trial_bits = 64
        trial = [trial_integer(proposal[i, 0], trial_bits) for i in range(count)]
        x = [arb(fmpq(i, 1 << trial_bits)) for i in trial]
        upper_debit = gram[0][0]-2*sum((x[i]*T[i] for i in range(count)), arb(0))
        upper_debit += sum((x[i]*K[i][j]*x[j] for i in range(count)
                            for j in range(count)), arb(0))
        margin = retained-upper_debit
        cases.append({'selected': list(selected[:count]), 'trial_integers': [str(i) for i in trial],
                      'trial_grid_bits': trial_bits, 'debit_interval': grid_interval(upper_debit, args.output_bits),
                      'margin_interval': grid_interval(margin, args.output_bits),
                      'positive_margin': bool(margin > 0),
                      'improvement_interval': grid_interval(gram[0][0]-upper_debit, args.output_bits)})
        print('exterior columns', count, 'debit', upper_debit.str(24),
              'margin', margin.str(24), flush=True)

    dependencies = ['arb_arithmetic.py', 'arb_far_coupling.py', 'arb_tail_moments.py',
                    'arb_middle_scalars.py', 'diagonal_floor_witness.py']
    result = {
        'schema': 'ACTUAL_C27_INVERSE_RESIDUAL_PILOT_V1', 'c': 27, 'N': 1280, 'M': far['M'],
        'basis': scalar['basis'], 'coefficients': [{'n': n, 'integer': value} for n, value in v.items()],
        'selected_exterior_columns': list(selected),
        'working_bits': args.precision, 'grid_bits': args.output_bits,
        'retained_quadratic': grid_interval(retained, args.output_bits),
        'old_full_debit': grid_interval(gram[0][0], args.output_bits),
        'joint_gram': [[grid_interval(x, args.output_bits) for x in row] for row in gram],
        'F': [[grid_interval(x, args.output_bits) for x in row] for row in F],
        'B': [grid_interval(x, args.output_bits) for x in B],
        'middle_weight_bits': weight_bits, 'middle_weight_numerators': weights,
        'far_weight_numerator': str(D.numerator), 'far_weight_denominator': str(D.denominator),
        'extended_feature_last_mode': max(selected), 'reused_moment_intervals': 168,
        'extended_far_matrix_error': grid_interval(far_error, args.output_bits),
        'cases': cases,
        'scalar_packet_sha256': SCALAR_SHA256, 'witness_sha256': WITNESS_SHA256,
        'middle_packet_sha256': MIDDLE_SHA256, 'far_packet_sha256': FAR_SHA256,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'dependency_sha256': {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                              for name in dependencies},
        'python_flint': flint.__version__, 'FLINT': flint.__FLINT_VERSION__,
        'elapsed_seconds': time.monotonic()-started,
        'conditions': [
            'The existing actual Fourier/form identification is assumed.',
            'The selected exterior Fourier basis vectors belong to D(H_ext).',
            'The same positive actual exterior operator satisfies H_ext^(-1) <= W '
            'by the existing even Neumann-Gamma and prime comparison.'],
        'scope': 'Directed efficacy test of the existing finite-column inverse-residual '
                 'comparison on the previously known high-frequency witness direction. '
                 'All middle and infinite far rows use the same actual columns and weights. '
                 'Positivity of the actual Schur quadratic on this vector already follows '
                 'from G7; this tests a sufficient estimate. No full1281-mode retained '
                 'sign, new general theorem, Lean or RH certificate.'}
    raw = (json.dumps(result, sort_keys=True, separators=(',', ':'))+'\n').encode()
    args.out.write_bytes(zlib.compress(raw, level=9))
    print('result_sha256', hashlib.sha256(args.out.read_bytes()).hexdigest(), flush=True)


if __name__ == '__main__':
    main()
