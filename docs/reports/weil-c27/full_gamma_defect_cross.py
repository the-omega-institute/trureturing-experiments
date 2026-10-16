"""Enclose the common lifted cross row from the accepted scalar defect.

The identity r=theta*u+e transports the saved residual's channel map
to u without regenerating any channel values. Existing weighted Cauchy
bounds pay e and the far cross; no inverse or scientific supplier runs.
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


def consume(args):
    started = time.monotonic()
    projection, projection_hash = read_json(args.projection)
    screen, screen_hash = read_json(args.screen)
    common, common_hash = read_json(args.common)
    boundary, boundary_hash = read_json(args.boundary)
    base, base_hash = read_json(args.base)
    defect, defect_hash = read_json(args.defect)
    direction, direction_hash = read_json(args.direction)
    core, core_hash = read_json(args.core)
    if (projection['schema'], screen['schema'], common['schema'], boundary['schema'],
            defect['schema'], direction['schema'], screen['projection_sha256'],
            screen['common_column_sha256'], screen['boundary_sha256'], screen['base_target_sha256'],
            projection['common_column_sha256'], projection['boundary_sha256'],
            projection['base_target_sha256'], common['finite_defect_sha256'],
            common['direction_sha256'], defect['metadata'], defect['metadata']['direction_sha256'],
            direction['updated_target_sha256'], direction['selected_remaining_mode']) != (
            'ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1', 'ACTUAL_C27_FULL_GAMMA_CROSS_SCREEN_V1',
            'ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1', 'ACTUAL_C27_JOINT_GAMMA_ACTION_V1',
            'ACTUAL_C27_FINITE_DUAL_DEFECT_V1', 'ACTUAL_C27_JOINT_GAMMA_DIRECTION_V1',
            projection_hash, common_hash, boundary_hash, base_hash, common_hash, boundary_hash,
            base_hash, defect_hash, direction_hash, common['metadata'], direction_hash, base_hash, 0):
        raise ValueError('Require the same accepted scalar defect, exact direction and old-base thin data.')
    if (core['schema'], core_hash, core['target_packet_sha256'], core['remaining_modes']) != (
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1', projection['base_core_sha256'],
            base_hash, projection['remaining_modes']):
        raise ValueError('Require the matching old-base complement for the scalar metric budgets.')
    expected = projection['input_sha256']
    if any(x != expected for x in (screen['input_sha256'], common['metadata']['input_sha256'],
                                  boundary['input_sha256'], base['metadata']['input_sha256'])):
        raise ValueError('Different scientific input provenance.')
    for name, digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Scientific input hash mismatch: '+name)
    modes = projection['remaining_modes']
    if screen['remaining_modes'] != modes or [r['mode'] for r in screen['rows']] != modes:
        raise ValueError('Different lift coordinate order.')
    c = matrix_in(projection['common_column_on_lift'], projection['grid_bits'])
    betaV = matrix_in(projection['channel_map_on_lift'], projection['grid_bits'])
    beta = matrix_in(boundary['finite_beta'], boundary['grid_bits'])
    K = matrix_in(base['kernel'], base['grid_bits'])
    rank, channels = len(projection['far_features_on_lift']), len(boundary['bins'])
    values = direction['direction_dyadic_integers']
    if (c.nrows(), c.ncols(), betaV.nrows(), betaV.ncols(), beta.nrows(), beta.ncols(),
            K.nrows(), K.ncols(), len(values), channels) != (
            1, len(modes), channels, len(modes), channels, N+1,
            rank+channels, rank+channels, N+1, 19):
        raise ValueError('Incomplete exact direction or thin channel maps.')
    v0 = arb_mat([[arb(fmpq(int(x), 1 << direction['direction_dyadic_bits']))] for x in values])
    theta = arb(fmpq(int(defect['theta_integer']), 1 << defect['theta_bits']))
    energy = interval_ball(common['bounds']['middle_action_energy'], common['grid_bits'])
    error = interval_ball(defect['weighted_residual_defect'], defect['grid_bits'])
    if not abs(theta) > 0 or not energy > 0 or not error >= 0:
        raise ArithmeticError('Require nonzero exact theta, positive action energy and nonnegative defect.')
    inverse = arb_mat([[-K[rank+i, rank+j] for j in range(channels)]
                       for i in range(channels)])
    # Fphi=(Beta0*v0)^*/theta-e^*W_J U_J/theta.
    center = c-(beta*v0).transpose()*inverse*betaV/theta
    factor = (error.upper()/(theta**2*energy.lower())).sqrt().upper()
    rho = arb(fmpq(1, 32))
    denominator = (interval_ball(common['bounds']['rank_one_delta_upper'], common['grid_bits']).upper()
                   +rho*interval_ball(common['bounds']['far_action_energy_upper'], common['grid_bits']).upper()).upper()
    schur = matrix_in(core['reduced_schur_intervals'], core['grid_bits'])
    rows, enclosure = [], []
    for j, old_row in enumerate(screen['rows']):
        finite_error = (interval_ball(old_row['finite_channel_cauchy_upper'],
                                     screen['grid_bits']).upper()*factor).upper()
        far_error = interval_ball(old_row['far_cross_cauchy_upper'], screen['grid_bits']).upper()
        radius = (finite_error+far_error).upper()
        value = center[0, j]+arb(0, radius)
        if not schur[j, j] < 0 or not denominator > 0:
            raise ArithmeticError('This scalar budget uses the old negative diagonal and positive charge.')
        squared_lower = abs(value).lower()**2
        gain_diagonal_lower = (schur[j, j].lower()+squared_lower/denominator).lower()
        metric_budget_lower = (squared_lower/(-schur[j, j].lower())-denominator).lower()
        enclosure.append(value)
        rows.append(dict(mode=modes[j],
                         finite_center=envelope(center[0, j], args.output_bits),
                         finite_defect_error_upper=envelope(finite_error, args.output_bits),
                         far_cross_error_upper=envelope(far_error, args.output_bits),
                         full_common_row=envelope(value, args.output_bits),
                         old_lift_gain_diagonal_lower=envelope(gain_diagonal_lower, args.output_bits),
                         scalar_good_metric_budget_lower=envelope(metric_budget_lower, args.output_bits),
                         nonzero_certified=bool(abs(value) > 0)))
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_DEFECT_CROSS_V1',
                  projection_sha256=projection_hash, screen_sha256=screen_hash,
                  common_column_sha256=common_hash, boundary_sha256=boundary_hash,
                  base_target_sha256=base_hash, scalar_defect_sha256=defect_hash,
                  exact_direction_sha256=direction_hash, input_sha256=expected,
                  base_core_sha256=core_hash,
                  remaining_modes=modes, rho='1/32', grid_bits=args.output_bits,
                  working_bits=args.precision, rows=rows,
                  common_row_on_lift=matrix_out(arb_mat([enclosure]), args.output_bits),
                  finite_error_contraction_factor=envelope(factor, args.output_bits),
                  common_denominator_upper=envelope(denominator, args.output_bits),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started,
                  conditions=screen['conditions']+[
                      'The saved exact scalar direction has r_J=R_J v0 and defect e=r_J-theta*(Mw)_J.',
                      'The accepted scalar defect encloses e^*W_J e for this same exact theta and original weights.'],
                  scope='Enclosures of the full common row on the old exact18-column lift only. '
                    'The scalar residual defect bounds the missing finite channel contraction; '
                    'the inherited full far cross bound remains paid. Common denominator '
                    'delta_up+rho*E_M is an upper for the inverse-residual square-completion charge '
                    'under the matching inverse comparison. Lift diagonal bounds pay this charge '
                    'but omit the good-block inverse metric. A scalar metric below its stated '
                    'budget would repair only that updated Schur diagonal, not the full matrix. '
                    'The metric is not supplied; no full coupled sign, new generic method, Lean or RH claim.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print('defect-transported common row nonzero', sum(x['nonzero_certified'] for x in rows),
          'factor', factor.str(16), 'seconds', time.monotonic()-started, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(projection='full-gamma-lift-projection.json.zlib',
                             screen='full-gamma-cross-screen.json',
                             common='full-gamma-common-column.json.zlib',
                             boundary='joint-gamma-6765-finite.json.zlib',
                             base='joint-gamma-6765-rho32-v2-target.json.zlib',
                             defect='full-gamma-dual-32768-defect.json.zlib',
                             direction='joint-gamma-6765-direction.json',
                             core='joint-gamma-6765-rho32-v2-schur.json').items():
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
