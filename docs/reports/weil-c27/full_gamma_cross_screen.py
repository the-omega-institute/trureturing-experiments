"""Screen the missing common cross row with cached weighted Cauchy bounds.

This consumes the thin lift packet, not a new target or an inverse action.
The finite channel contraction remains unknown. An inconclusive triangle
bound is not a feasible cancellation or an obstruction to a stronger bound.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from flint import arb, arb_mat, ctx, fmpq

from arb_inverse_residual_obstruction import HERE, M, N
from diagonal_floor_witness import interval_ball
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in
from structured_inverse_residual_interval import envelope


def consume(args):
    started = time.monotonic()
    projection, projection_hash = read_json(args.projection)
    common, common_hash = read_json(args.common)
    boundary, boundary_hash = read_json(args.boundary)
    base, base_hash = read_json(args.base)
    if (projection['schema'], common['schema'], boundary['schema'], base['schema'],
            projection['common_column_sha256'], projection['boundary_sha256'],
            projection['base_target_sha256'], common['target_sha256'],
            boundary['target_packet_sha256'], base['metadata']['original_target_sha256'],
            base['metadata']['joint_boundary_action_sha256'],
            base['metadata']['common_young_parameter'], projection['N'], projection['M']) != (
            'ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
            'ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1',
            'ACTUAL_C27_JOINT_GAMMA_ACTION_V1', 'ACTUAL_C27_STRUCTURED_TARGET_V1',
            common_hash, boundary_hash, base_hash, projection['original_target_sha256'],
            projection['original_target_sha256'], projection['original_target_sha256'],
            boundary_hash, '1/32', N, M):
        raise ValueError('Require the accepted same-source common column and old19-channel lift.')
    expected = projection['input_sha256']
    if any(x != expected for x in (common['metadata']['input_sha256'],
                                  boundary['input_sha256'], base['metadata']['input_sha256'])):
        raise ValueError('Different scientific input provenance.')
    for name, digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Scientific source hash mismatch: '+name)
    c = matrix_in(projection['common_column_on_lift'], projection['grid_bits'])
    beta = matrix_in(projection['channel_map_on_lift'], projection['grid_bits'])
    far = matrix_in(projection['far_gram_upper'], projection['grid_bits'])
    A0 = matrix_in(boundary['finite_gram'], boundary['grid_bits'])
    K = matrix_in(base['kernel'], base['grid_bits'])
    rank = len(projection['far_features_on_lift'])
    channels, count = len(boundary['bins']), len(projection['remaining_modes'])
    if (c.nrows(), c.ncols(), beta.nrows(), beta.ncols(), far.nrows(), far.ncols(),
            A0.nrows(), A0.ncols(), K.nrows(), K.ncols(), channels) != (
            1, count, channels, count, count, count, channels, channels,
            rank+channels, rank+channels, 19):
        raise ValueError('Incomplete old-base inverse kernel or projected common data.')
    # The saved base already contains -A_rho^-1 as its19-channel kernel.
    # This is not the positive good block's inverse or a fresh solve.
    inverse = arb_mat([[-K[rank+i, rank+j] for j in range(channels)]
                       for i in range(channels)])
    Z = inverse*beta
    energies = Z.transpose()*A0*Z
    middle_energy = interval_ball(common['bounds']['middle_action_energy'], common['grid_bits'])
    far_energy = interval_ball(common['bounds']['far_action_energy_upper'], common['grid_bits'])
    if not middle_energy >= 0 or not far_energy >= 0:
        raise ArithmeticError('Nonnegative effective-action energy bounds are required.')
    rho = arb(fmpq(1, 32))
    rows = []
    for j, mode in enumerate(projection['remaining_modes']):
        if not energies[j, j] > 0 or not far[j, j] > 0:
            raise ArithmeticError('This screen requires positive enclosed squared bounds.')
        finite_radius = (middle_energy.upper()*energies[j, j].upper()).sqrt().upper()
        far_radius = ((1+rho)*(far_energy.upper()*far[j, j].upper()).sqrt()).upper()
        c_abs_lower = abs(c[0, j]).lower()
        raw_margin = c_abs_lower-far_radius
        complete_margin = raw_margin-finite_radius
        rows.append(dict(mode=mode,
                         raw_column_abs_lower=envelope(c_abs_lower, args.output_bits),
                         finite_channel_cauchy_upper=envelope(finite_radius, args.output_bits),
                         far_cross_cauchy_upper=envelope(far_radius, args.output_bits),
                         raw_column_minus_far_lower=envelope(raw_margin.lower(), args.output_bits),
                         complete_triangle_lower=envelope(complete_margin.lower(), args.output_bits),
                         raw_column_survives_far=bool(raw_margin > 0),
                         complete_row_nonzero_certified=bool(complete_margin > 0)))
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_CROSS_SCREEN_V1',
                  projection_sha256=projection_hash, common_column_sha256=common_hash,
                  boundary_sha256=boundary_hash, base_target_sha256=base_hash,
                  input_sha256=expected, remaining_modes=projection['remaining_modes'],
                  grid_bits=args.output_bits, working_bits=args.precision, rho='1/32',
                  rows=rows, elapsed_seconds=time.monotonic()-started,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  conditions=projection['conditions']+[
                      'A0 is the Gram of all middle channel rows, and hence majorizes the J-only channel Gram.',
                      'The old base channel kernel is exactly -A_rho^-1 for the same Beta0.',
                      'The common far cross uses the same Mw and R and the weight (1+rho)W_far.'],
                  scope='Per-coordinate weighted Cauchy and triangle bounds only. '
                    'A positive raw margin excludes cancellation of b_fin^*V by the far cross alone. '
                    'Only a positive complete margin certifies a nonzero coordinate of '
                    'b_fin^*V-Fphi A_rho^-1 BetaV+(1+rho)(Mw)_far^*W_far R_far V. '
                    'An inconclusive complete margin neither constructs feasible joint cancellation '
                    'nor excludes sharper estimates. No common monotone gain, good-block inverse '
                    'metric, matrix sign, Lean verification, novelty or RH result is supplied.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print('raw coordinates surviving far', sum(x['raw_column_survives_far'] for x in rows),
          'complete coordinates certified', sum(x['complete_row_nonzero_certified'] for x in rows),
          'seconds', time.monotonic()-started, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(projection='full-gamma-lift-projection.json.zlib',
                             common='full-gamma-common-column.json.zlib',
                             boundary='joint-gamma-6765-finite.json.zlib',
                             base='joint-gamma-6765-rho32-v2-target.json.zlib').items():
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
