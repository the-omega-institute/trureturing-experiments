"""Transport new common data through the saved nineteen-channel Schur lift.

Only thin products are missing. Reuse the exact implicit lift enclosed by
the validated inverse action, the new common column, old channel map and
signed far majorant. No target reconstruction, inverse solve or supplier.
"""
import argparse
import hashlib
from pathlib import Path
import time

from flint import arb, arb_mat, ctx, fmpq

from arb_inverse_residual_obstruction import HERE, M, N
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, matrix_out, original_target, write_packet
from structured_residual_schur_finish import load_target


def consume(args):
    started = time.monotonic()
    old, old_hash, metadata_bits = original_target(args)
    base, prefix, remaining, base_hash, pivot_raw = load_target(args.lift_target, args.lift_pivots)
    actions, actions_hash = read_json(args.actions)
    core, core_hash = read_json(args.core)
    boundary, boundary_hash = read_json(args.boundary)
    common, common_hash = read_json(args.common)
    pivot_hash = hashlib.sha256(pivot_raw).hexdigest()
    if (base['metadata']['original_target_sha256'],
            base['metadata']['joint_boundary_action_sha256'],
            base['metadata']['common_young_parameter'],
            boundary['schema'], boundary['target_packet_sha256'],
            common['schema'], common['target_sha256']) != (
            old_hash, boundary_hash, '1/32', 'ACTUAL_C27_JOINT_GAMMA_ACTION_V1',
            old_hash, 'ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1', old_hash):
        raise ValueError('Require the original residual, accepted19-channel base and matching new column.')
    if (actions['schema'], actions['target_packet_sha256'], actions['pivot_result_sha256'],
            actions['positive_modes'], actions['remaining_modes'], core['schema'],
            core['target_packet_sha256'], core['pivot_result_sha256'], core['remaining_modes'],
            core['reused_positive_coordinate_count']) != (
            'ACTUAL_C27_SCHUR_LIFT_V1', base_hash, pivot_hash, prefix, remaining,
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1', base_hash, pivot_hash, remaining, len(prefix)):
        raise ValueError('Require the complete exact implicit lift and core of this same base.')
    expected = old['metadata']['input_sha256']
    if any(x != expected for x in (base['metadata']['input_sha256'], boundary['input_sha256'],
                                   common['metadata']['input_sha256'], core['input_sha256'])):
        raise ValueError('Different actual scientific source inputs.')
    if (base['metadata']['original_metadata_grid_bits'], metadata_bits,
            common['metadata']['N'], common['metadata']['M'],
            common['metadata']['last_mode'], common['metadata']['basis']) != (
            metadata_bits, 512, N, M, 2048, old['metadata']['basis']):
        raise ValueError('Require inherited512-bit allowances and the same complete residual basis.')
    lift = matrix_in(actions['inverse_action'], actions['grid_bits'])
    column = matrix_in(common['column'], common['grid_bits'])
    beta = matrix_in(boundary['finite_beta'], boundary['grid_bits'])
    G = matrix_in(old['generators'], old['grid_bits'])
    K = matrix_in(old['kernel'], old['grid_bits'])
    count, channels, rank = len(remaining), len(boundary['bins']), K.nrows()
    if (lift.nrows(), lift.ncols(), column.nrows(), column.ncols(), beta.nrows(), beta.ncols(),
            G.nrows(), G.ncols(), K.ncols(), channels) != (
            len(prefix), count, N+1, 1, channels, N+1, N+1, 4+2*rank, rank, 19):
        raise ValueError('Incomplete saved lift, column, channel map or signed far factors.')
    # V=[-P^-1 B_g;I] in the saved coordinate order. The enclosed exact
    # inverse action is used directly, without rounding to a dyadic lift.
    V = arb_mat(N+1, count)
    for i, n in enumerate(prefix):
        for j in range(count):
            V[n, j] = -lift[i, j]
    for j, n in enumerate(remaining):
        V[n, j] = 1
    P = arb_mat([[G[n, 4+j] for j in range(rank)] for n in range(N+1)])
    projected_features = P.transpose()*V
    lift_gram = V.transpose()*V
    # Reuse the upward allowance already transported by this base. Its
    # explicitly labelled grid is distinct from the original metadata grid.
    eps_bits = base['metadata']['reused_far_allowance_grid_bits']
    eps = fmpq(int(base['metadata']['reused_far_allowance'][1]), 1 << eps_bits)
    original_eps = fmpq(int(old['metadata']['complete_far_operator_allowance'][1]),
                        1 << metadata_bits)
    if eps < original_eps:
        raise ArithmeticError('The transported far majorant must retain the original full allowance.')
    far_gram = projected_features.transpose()*K*projected_features+arb(eps)*lift_gram
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
                  N=N, M=M, last_mode=2048, basis=old['metadata']['basis'],
                  original_target_sha256=old_hash, base_target_sha256=base_hash,
                  base_pivots_sha256=pivot_hash, base_core_sha256=core_hash,
                  base_lift_sha256=actions_hash, boundary_sha256=boundary_hash,
                  common_column_sha256=common_hash, input_sha256=expected,
                  positive_modes=prefix, remaining_modes=remaining,
                  original_metadata_grid_bits=metadata_bits,
                  exact_far_allowance=str(eps), grid_bits=args.output_bits,
                  working_bits=args.precision,
                  common_column_on_lift=matrix_out(column.transpose()*V, args.output_bits),
                  channel_map_on_lift=matrix_out(beta*V, args.output_bits),
                  far_features_on_lift=matrix_out(projected_features, args.output_bits),
                  lift_gram=matrix_out(lift_gram, args.output_bits),
                  far_gram_upper=matrix_out(far_gram, args.output_bits),
                  effective_action_far_energy_upper=common['bounds']['far_action_energy_upper'],
                  effective_action_far_energy_grid_bits=common['grid_bits'],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started,
                  conditions=common['conditions']+[
                      'The saved19-channel base lift encloses its exact implicit P^-1 B_g.',
                      'The original signed far factors plus their full allowance majorize the same R_far^* W R_far.'],
                  scope='Only shared thin data transported to the existing exact18-column lift. '
                    'The signed core and full allowance are retained. No new inverse, scalar or moment supplier, '
                    'special-function/channel production or target sign test. The actual common far cross row, '
                    'monotone gain and good-block inverse metric remain missing. No Lean, novelty or RH claim.')
    write_packet(args.out, result)
    print('saved lift projections enclosed', count, 'channels', channels,
          'seconds', time.monotonic()-started, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    defaults = dict(target='structured-joint-2048-target.json.zlib',
                    pivots='structured-joint-2048-pivots.json',
                    lift_target='joint-gamma-6765-rho32-v2-target.json.zlib',
                    lift_pivots='joint-gamma-6765-rho32-v2-pivots.json',
                    actions='joint-gamma-6765-rho32-v2-actions.json.zlib',
                    core='joint-gamma-6765-rho32-v2-schur.json',
                    boundary='joint-gamma-6765-finite.json.zlib',
                    common='full-gamma-common-column.json.zlib')
    for name, default in defaults.items():
        parser.add_argument('--'+name.replace('_', '-'), type=Path, default=HERE/default)
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.output_bits != 640 or args.precision < 832:
        parser.error('Use a fresh output, inherited640-bit grid and at least192 guard bits.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
