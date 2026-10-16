"""Calibrate a proposed source quadratic before any new source search.

Use cached two-source rows only as interface inputs. The fixed-linear
variational rule differs from the archived fixed rational trials. Exact
rational arithmetic checks the two expressions on one midpoint model;
ordinary floating evaluation and interval contraction are compared with
that same model. No inverse action or scientific supplier is regenerated.
"""
import argparse
import hashlib
import platform
from pathlib import Path
import time

import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits
import flint
from flint import arb, ctx, fmpq, fmpq_mat

from arb_inverse_residual_obstruction import HERE
from full_gamma_correlated_sign import symmetric
from full_gamma_directional_metric import balls, midpoint
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope


def point(rows, bits):
    return fmpq_mat([[midpoint(pair, bits) for pair in row] for row in rows])


def sym(matrix):
    return (matrix+matrix.transpose())/2


def coefficient(sigma, gram, reference, kappa):
    """The known variational/square-completion formula, in either field."""
    n = sigma.nrows()
    kind = type(sigma)
    a = sigma[0, 0]
    e = kind([[int(i == 0)] for i in range(n)])
    T = kind([[-sigma[0, j+1]/a if i == 0 else int(i == j+1)
               for j in range(n-1)] for i in range(n)])
    P = sym(T.transpose()*sigma*T)
    A = kappa-kappa**2*gram[0, 0]/(-a)
    L = e*e.transpose()/a+T*reference*T.transpose()
    Q = (e*e.transpose()/a+T*(2*reference-reference*P*reference)*T.transpose()
         -kappa**2/A*L.transpose()*gram*L)
    return sym(Q), T, P, A, a, e


def factored(row, charge, T, P, A, a, e, gram, reference, kappa, trial=None):
    H = row*T
    X = reference*H.transpose() if trial is None else trial
    b0 = row*e
    g = e*b0.transpose()/a+T*X
    C = (charge+b0*b0.transpose()/a+H*X+X.transpose()*H.transpose()
         -X.transpose()*P*X-kappa**2/A*g.transpose()*gram*g)
    return sym(C), X


def floats(matrix):
    return np.array([[float(matrix[i, j])
                      for j in range(matrix.ncols())] for i in range(matrix.nrows())])


def float_exact(value):
    numerator, denominator = float(value).as_integer_ratio()
    return fmpq(numerator, denominator)


def produce(args):
    started = time.monotonic()
    paths = dict(two=HERE/'full-gamma-two-source.json.zlib',
                 core=HERE/'joint-gamma-6765-rho32-v2-schur.json',
                 projection=HERE/'full-gamma-lift-projection.json.zlib',
                 reference=HERE/'full-gamma-directional-metric.json.zlib')
    packets, hashes = {}, {}
    for name, path in paths.items():
        packets[name], hashes[name] = read_json(path)
    two, core, projection, reference = (packets[name] for name in paths)
    if (two['schema'], core['schema'], projection['schema'], reference['schema']) != (
            'ACTUAL_C27_TWO_SOURCE_JOINT_LOWER_V1',
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1',
            'ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
            'ACTUAL_C27_FULL_GAMMA_DIRECTIONAL_METRIC_V1'):
        raise ValueError('Require the accepted same-source cached model.')
    if any(two['packet_sha256'][name] != hashes[name] for name in ('core', 'projection')):
        raise ValueError('Different source/lift provenance.')
    if two['packet_sha256']['certificate'] != hashes['reference']:
        raise ValueError('Different exact variational reference.')
    kq = fmpq(33, 32)
    eps = fmpq(projection['exact_far_allowance'])
    sq = sym(point(core['reduced_schur_intervals'], core['grid_bits'])
             +2*kq*eps*point(projection['lift_gram'], projection['grid_bits']))
    gq = sym(point(projection['far_gram_upper'], projection['grid_bits']))
    rq = fmpq_mat([[fmpq(x) for x in row] for row in reference['reference_inverse']])
    bq = point(two['finite_row_on_lift'], two['grid_bits'])
    kjq = sym(point(two['finite_charge_matrix'], two['grid_bits']))
    if rq != rq.transpose() or bq.nrows() != 2 or bq.ncols() != 18:
        raise ValueError('Require a symmetric exact reference and the complete old two-source rows.')
    Qq, Tq, Pq, Aq, aq, eq = coefficient(sq, gq, rq, kq)
    if not Aq > 0 or not aq < 0:
        raise ArithmeticError('Different negative-pivot/far model.')
    exact, Xq = factored(bq, kjq, Tq, Pq, Aq, aq, eq, gq, rq, kq)
    if exact != sym(kjq+bq*Qq*bq.transpose()):
        raise ArithmeticError('The exact rational quadratic identity failed.')

    with threadpool_limits(limits=1):
        B, K, Q = floats(bq), floats(kjq), floats(Qq)
        T, P, R, G, e = map(floats, (Tq, Pq, rq, gq, eq))
        a, A, kappa = map(float, (aq, Aq, kq))
        expanded = K+B@Q@B.T
        H, b0 = B@T, B@e
        X = R@H.T
        grad = e@b0.T/a+T@X
        linear = K+b0@b0.T/a+H@X+X.T@H.T-X.T@P@X-kappa**2/A*grad.T@G@grad
        frozen_X = floats(Xq)
        frozen_grad = e@b0.T/a+T@frozen_X
        frozen = (K+b0@b0.T/a+H@frozen_X+frozen_X.T@H.T
                  -frozen_X.T@P@frozen_X-kappa**2/A*frozen_grad.T@G@frozen_grad)
        runtime = [dict(user_api=info['user_api'], internal_api=info['internal_api'],
                        version=info['version'], num_threads=info['num_threads'])
                   for info in threadpool_info()]
    variants = dict(expanded_best_case_rounded_exact_coefficient=expanded,
                    factorized_float_reference=linear,
                    factorized_rounded_exact_trial=frozen)
    if any(not np.isfinite(value).all() for value in variants.values()):
        raise ArithmeticError('Nonfinite calibration; no source selection is permitted.')
    errors = {name: max(abs(float_exact(value[i, j])-exact[i, j])
                        for i in range(2) for j in range(2))
              for name, value in variants.items()}
    scale = max(fmpq(1), *(abs(x) for x in exact.entries()))
    tolerance = scale*fmpq(1, 100000000)

    sb = symmetric(matrix_in(core['reduced_schur_intervals'], core['grid_bits'])
                   +2*arb(kq)*arb(eps)*matrix_in(projection['lift_gram'], projection['grid_bits']))
    gb = symmetric(matrix_in(projection['far_gram_upper'], projection['grid_bits']))
    bb = matrix_in(two['finite_row_on_lift'], two['grid_bits'])
    kjb = symmetric(matrix_in(two['finite_charge_matrix'], two['grid_bits']))
    Qb, Tb, Pb, Ab, ab, eb = coefficient(sb, gb, balls(rq), arb(kq))
    interval_factored, _ = factored(bb, kjb, Tb, Pb, Ab, ab, eb, gb, balls(rq), arb(kq))
    interval_frozen, _ = factored(bb, kjb, Tb, Pb, Ab, ab, eb, gb, balls(rq), arb(kq),
                                 trial=balls(Xq))
    interval_expanded = symmetric(kjb+bb*Qb*bb.transpose())
    for matrix in (interval_factored, interval_frozen, interval_expanded):
        if any(not matrix[i, j].contains(arb(exact[i, j])) for i in range(2) for j in range(2)):
            raise ArithmeticError('A reference point is outside its interval model.')
    width = lambda matrix: max(x.upper()-x.lower() for x in matrix.entries())
    scalars = dict(reference_scale=arb(scale), calibration_tolerance=arb(tolerance),
                   factorized_max_width=width(interval_factored),
                   frozen_trial_max_width=width(interval_frozen),
                   expanded_max_width=width(interval_expanded),
                   exact_trial_max_abs=max(abs(arb(x)) for x in Xq.entries()),
                   exact_coefficient_max_abs=max(abs(arb(x)) for x in Qq.entries()))
    scalars.update({name+'_max_abs_error': arb(error) for name, error in errors.items()})
    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result = dict(schema='ACTUAL_C27_SOURCE_QUADRATIC_CALIBRATION_V1',
                  packet_sha256=hashes, source_sha256=source_hash,
                  grid_bits=args.output_bits, working_bits=args.precision,
                  exact_rational_identity=True, exact_midpoint_trial_rule='Rref*Tstar*b',
                  frozen_interval_trial_rule='The same reference Xq is fixed exactly; it is not the variable linear-trial operator.',
                  runtime=dict(python=platform.python_version(), flint=flint.__version__,
                               numpy=np.__version__, machine=platform.machine(), libraries=runtime),
                  archived_trial_recomputed=False, new_inverse_actions=0,
                  new_scalar_or_channel_rows=0, eigen_operator_actions=0,
                  float_pass={name: error <= tolerance for name, error in errors.items()},
                  float_values={name: value.tolist() for name, value in variants.items()},
                  reference_matrix=matrix_out(balls(exact), args.output_bits),
                  interval_factorized=matrix_out(interval_factored, args.output_bits),
                  interval_frozen_trial=matrix_out(interval_frozen, args.output_bits),
                  interval_expanded=matrix_out(interval_expanded, args.output_bits),
                  scalar_bounds={name: envelope(value, args.output_bits) for name, value in scalars.items()},
                  verdict=('ordinary-double-source-proposal-calibration-failed'
                           if any(error > tolerance for name, error in errors.items()
                                  if name != 'factorized_rounded_exact_trial') else 'calibration-passed'),
                  elapsed_seconds=time.monotonic()-started,
                  scope='Necessary implementation calibration only. The exact reference is the stated '
                    'rational midpoint fixed-linear-trial model, not the archived fixed-trial C. '
                    'A failed gate stops this ordinary-double proposal, not mathematical source existence. '
                    'Interval widths measure this evaluation form, not intrinsic input uncertainty. '
                    'No new source, positivity/exclusion theorem, generic algorithm, Lean or RH result.')
    for name, path in paths.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != hashes[name]:
            raise ValueError('A consumed cache changed: '+name)
    write_packet(args.out, result)
    for name, value in scalars.items():
        print(name, value.str(14), flush=True)
    print(result['verdict'], result['float_pass'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib' or args.output_bits != 640 or args.precision < 832:
        parser.error('Use fresh .zlib output and the inherited640-bit grid with192 guardbits.')
    ctx.prec = args.precision
    produce(args)


if __name__ == '__main__':
    main()
