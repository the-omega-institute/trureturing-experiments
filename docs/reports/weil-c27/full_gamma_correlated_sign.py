"""Bound the full lifted rank-one sign using the shared error Grams.

Reuse the old one-negative core and scalar-defect transport. A standard
validated 17-dimensional solve contracts the positive complement; no
scientific supplier, large inverse, determinant or eigensolver runs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from flint import arb, arb_mat, ctx, fmpq

from arb_inverse_residual_obstruction import HERE, N
from diagonal_floor_witness import interval_ball
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, matrix_out
from structured_inverse_residual_interval import envelope


def symmetric(matrix):
    """Mirror an enclosing upper triangle of the same exact symmetric matrix."""
    if matrix.nrows() != matrix.ncols():
        raise ValueError('A square real symmetric enclosure is required.')
    n = matrix.nrows()
    for i in range(n):
        for j in range(i):
            if not matrix[i, j].overlaps(matrix[j, i]):
                raise ArithmeticError('Incompatible symmetric entry enclosures.')
    return arb_mat([[matrix[min(i, j), max(i, j)] for j in range(n)] for i in range(n)])


def consume(args):
    started = time.monotonic()
    packets, hashes = {}, {}
    for name in ('projection', 'cross', 'common', 'boundary', 'base', 'defect',
                 'direction', 'core', 'rank'):
        packets[name], hashes[name] = read_json(getattr(args, name))
    projection, cross, common, boundary, base, defect, direction, core, rank = (
        packets[name] for name in ('projection', 'cross', 'common', 'boundary',
                                  'base', 'defect', 'direction', 'core', 'rank'))
    schemas = dict(projection='ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
                   cross='ACTUAL_C27_FULL_GAMMA_DEFECT_CROSS_V1',
                   common='ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1',
                   boundary='ACTUAL_C27_JOINT_GAMMA_ACTION_V1',
                   defect='ACTUAL_C27_FINITE_DUAL_DEFECT_V1',
                   direction='ACTUAL_C27_JOINT_GAMMA_DIRECTION_V1',
                   core='ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1',
                   rank='ACTUAL_C27_SMALL_SCHUR_RANK_FLOOR_V1')
    if any(packets[name]['schema'] != schema for name, schema in schemas.items()):
        raise ValueError('Require accepted old-base, common-row and inertia schemas.')
    for name in ('projection', 'cross'):
        if any(packets[name][key] != hashes[source] for key, source in (
                ('common_column_sha256', 'common'), ('boundary_sha256', 'boundary'),
                ('base_target_sha256', 'base'), ('base_core_sha256', 'core'))):
            raise ValueError('Different common row, boundary, old base or exact lift core.')
    if (cross['projection_sha256'], cross['scalar_defect_sha256'],
            cross['exact_direction_sha256'], common['finite_defect_sha256'],
            common['direction_sha256'], defect['metadata'],
            defect['metadata']['direction_sha256'], direction['updated_target_sha256'],
            core['target_packet_sha256'], rank['target_packet_sha256'],
            rank['schur_result_sha256'], rank['pivot_result_sha256'],
            projection['base_pivots_sha256'], base['metadata']['common_young_parameter']) != (
            hashes['projection'], hashes['defect'], hashes['direction'], hashes['defect'],
            hashes['direction'], common['metadata'], hashes['direction'], hashes['base'],
            hashes['base'], hashes['base'], hashes['core'], core['pivot_result_sha256'],
            core['pivot_result_sha256'], '1/32'):
        raise ValueError('Inconsistent scalar defect, exact direction or inherited inertia evidence.')
    expected = projection['input_sha256']
    if any(x != expected for x in (cross['input_sha256'], common['metadata']['input_sha256'],
                                  boundary['input_sha256'], base['metadata']['input_sha256'],
                                  core['input_sha256'])):
        raise ValueError('Different scientific input provenance.')
    for name, digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Scientific input hash mismatch: '+name)
    modes = projection['remaining_modes']
    if (len(modes), modes[0], rank['verdict'], rank['negative_count_lower'],
            rank['reused_positive_coordinate_count'], len(rank['leading_minors'])) != (
            18, 0, 'full-complement-inertia-counted', 1, 1263, 18):
        raise ValueError('Require the accepted old18-core with exactly one negative direction.')
    if any(packets[name]['remaining_modes'] != modes for name in ('cross', 'core', 'rank')):
        raise ValueError('Different old-lift coordinate order.')
    for k, minor in enumerate(rank['leading_minors'], 1):
        if (minor['order'], minor['determinant_sign'], minor['pivot_sign']) != (
                k, -1, -1 if k == 1 else 1) or int(minor['determinant_interval'][1]) >= 0:
            raise ValueError('The saved Jacobi signs do not establish the claimed positive complement.')
    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], core['grid_bits']))
    a = sigma[0, 0]
    if not a < 0:
        raise ArithmeticError('Require the accepted negative first pivot.')
    rest = len(modes)-1
    T = arb_mat([[(-sigma[0, j+1]/a) if i == 0 else arb(int(i == j+1))
                  for j in range(rest)] for i in range(rest+1)])
    positive = symmetric(arb_mat([[sigma[i+1, j+1]-sigma[i+1, 0]*sigma[0, j+1]/a
                                   for j in range(rest)] for i in range(rest)]))
    c = matrix_in(projection['common_column_on_lift'], projection['grid_bits'])
    betaV = matrix_in(projection['channel_map_on_lift'], projection['grid_bits'])
    beta = matrix_in(boundary['finite_beta'], boundary['grid_bits'])
    A0 = symmetric(matrix_in(boundary['finite_gram'], boundary['grid_bits']))
    kernel = matrix_in(base['kernel'], base['grid_bits'])
    features, channels = len(projection['far_features_on_lift']), len(boundary['bins'])
    values = direction['direction_dyadic_integers']
    if (c.nrows(), c.ncols(), betaV.nrows(), betaV.ncols(), beta.nrows(), beta.ncols(),
            A0.nrows(), kernel.nrows(), kernel.ncols(), len(values), channels) != (
            1, 18, channels, 18, channels, N+1, channels,
            features+channels, features+channels, N+1, 19):
        raise ValueError('Incomplete same-source thin maps or exact scalar direction.')
    inverse = symmetric(arb_mat([[-kernel[features+i, features+j] for j in range(channels)]
                                 for i in range(channels)]))
    v0 = arb_mat([[arb(fmpq(int(x), 1 << direction['direction_dyadic_bits']))] for x in values])
    theta = arb(fmpq(int(defect['theta_integer']), 1 << defect['theta_bits']))
    if not abs(theta) > 0:
        raise ArithmeticError('A nonzero exact scalar direction scale is required.')
    Ee = interval_ball(defect['weighted_residual_defect'], defect['grid_bits']).upper()
    EM = interval_ball(common['bounds']['far_action_energy_upper'], common['grid_bits']).upper()
    rho = arb(fmpq(1, 32))
    Z = inverse*betaV
    center = c-(beta*v0).transpose()*Z/theta
    Cf = symmetric((Ee/theta**2)*(Z.transpose()*A0*Z))
    Ct = symmetric((1+rho)**2*EM*matrix_in(projection['far_gram_upper'], projection['grid_bits']))
    h = center*T
    finite = symmetric(T.transpose()*Cf*T)
    far = symmetric(T.transpose()*Ct*T)
    # One standard validated solve with 35 thin right-hand sides.
    rhs = arb_mat([[h[0, i] if j == 0 else finite[i, j-1] if j <= rest
                    else far[i, j-rest-1] for j in range(1+2*rest)] for i in range(rest)])
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_CORRELATED_SIGN_V1',
                  packet_sha256=hashes, input_sha256=expected, remaining_modes=modes,
                  grid_bits=args.output_bits, working_bits=args.precision,
                  positive_complement_dimension=rest, solve_rhs_count=1+2*rest,
                  finite_error_gram_on_lift=matrix_out(Cf, args.output_bits),
                  far_error_gram_on_lift=matrix_out(Ct, args.output_bits),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  conditions=list(dict.fromkeys(cross['conditions']+rank['conditions']))+[
                      'The accepted old real symmetric core has inertia17positive,1negative,0zero.',
                      'Both error Grams bound components of the same actual full common row on the same exact lift.',
                      'The standard validated solve encloses the positive-complement inverse actions.'],
                  scope='Shared error-Gram contractions for a particular rank-one sufficient lower target. '
                    'Trace bounds and marginal scalar endpoints are conservative simultaneous bounds, '
                    'not claims of jointly attained optimizers. The missing good-block metric remains unpaid. '
                    'No actual Schur sign, cofinal positivity, new general method, Lean or RH claim.')
    try:
        actions = positive.solve(rhs, algorithm='precond')
        if any(not actions[i, j].is_finite() for i in range(rest) for j in range(1+2*rest)):
            raise ArithmeticError('The validated small solve returned a nonfinite enclosure.')
        center_sq = sum((h[0, i]*actions[i, 0] for i in range(rest)), arb(0))
        finite_sq = sum((actions[i, i+1] for i in range(rest)), arb(0))
        far_sq = sum((actions[i, i+rest+1] for i in range(rest)), arb(0))
        if any(not x.is_finite() or not x.upper() >= 0 for x in (center_sq, finite_sq, far_sq)):
            raise ArithmeticError('A nonnegative metric upper was not enclosed.')
        nc_lo, nc_up = max(arb(0), center_sq.lower()).sqrt().lower(), center_sq.upper().sqrt().upper()
        nf, nt = finite_sq.upper().sqrt().upper(), far_sq.upper().sqrt().upper()
        low_norm = max(arb(0), (nc_lo-nf-nt).lower())
        pos_lo, pos_up = (low_norm**2).lower(), ((nc_up+nf+nt)**2).upper()
        b0 = interval_ball(cross['common_row_on_lift'][0][0], cross['grid_bits'])
        negative = b0**2/a
        qlo, qup = (negative+pos_lo).lower(), (negative+pos_up).upper()
        D_interval = interval_ball(cross['common_denominator_upper'], cross['grid_bits'])
        Dlo, D = D_interval.lower(), D_interval.upper()
        if not Dlo > 0:
            raise ArithmeticError('Require a strictly positive common denominator upper.')
        slack = (-D-qup).lower()
        verdict = ('positive-good-metric-budget' if slack > 0 else
                   'Dstar-target-excluded-even-at-zero-good-metric' if qlo+Dlo > 0 else
                   'inconclusive-shared-trace-enclosure')
        scalars = dict(center_positive_energy=center_sq, finite_error_trace=finite_sq,
                       far_error_trace=far_sq, center_positive_norm_lower=nc_lo,
                       center_positive_norm_upper=nc_up, finite_error_norm_upper=nf,
                       far_error_norm_upper=nt, positive_energy_lower=pos_lo,
                       positive_energy_upper=pos_up, negative_energy=negative,
                       signed_inverse_metric_lower=qlo, signed_inverse_metric_upper=qup,
                       common_denominator_lower=Dlo, common_denominator_upper=D,
                       sufficient_good_metric_budget_lower=slack,
                       Dstar_target_scalar_lower=(qlo+Dlo).lower(),
                       Dstar_target_scalar_upper=(qup+D).upper())
        result.update(verdict=verdict, scalar_bounds={name: envelope(x, args.output_bits)
                                                    for name, x in scalars.items()})
        for name, value in scalars.items():
            print(name, value.str(14), flush=True)
    except (ArithmeticError, ValueError, ZeroDivisionError) as error:
        result.update(verdict='inconclusive-validated-contraction', diagnostic=str(error))
    result['elapsed_seconds'] = time.monotonic()-started
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(projection='full-gamma-lift-projection.json.zlib',
                             cross='full-gamma-defect-cross.json',
                             common='full-gamma-common-column.json.zlib',
                             boundary='joint-gamma-6765-finite.json.zlib',
                             base='joint-gamma-6765-rho32-v2-target.json.zlib',
                             defect='full-gamma-dual-32768-defect.json.zlib',
                             direction='joint-gamma-6765-direction.json',
                             core='joint-gamma-6765-rho32-v2-schur.json',
                             rank='joint-gamma-6765-rho32-v2-rank-floor.json').items():
        parser.add_argument('--'+name, type=Path, default=HERE/default)
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.output_bits != 640 or args.precision < 832:
        parser.error('Use a fresh output and the inherited640-bit grid with192 guard bits.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
