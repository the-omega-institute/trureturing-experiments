"""Consume the missing common retained column for a full-Gamma update.

Reuse the archived effective dual action and residual generators. The
existing finite Cauchy product computes R_J^*(W_J M_J w-w) on all retained
coordinates. No supplier, inverse, optimizer or retained quadratic runs.
"""
import argparse
import hashlib
from pathlib import Path
import time

from flint import arb, arb_mat, ctx, fmpq

from arb_inverse_residual_obstruction import HERE, M, N
from diagonal_floor_witness import interval_ball
from gamma_boundary_update import read_json
from gamma_channel_resolution import cauchy_projection
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope


def consume(args):
    started = time.monotonic()
    finite, finite_hash = read_json(args.finite)
    defect, defect_hash = read_json(args.defect)
    primal, primal_hash = read_json(args.primal)
    target, target_hash = read_json(args.target)
    generators, generator_hash = read_json(args.generators)
    pilot, pilot_hash = read_json(args.pilot)
    metadata = finite['metadata']
    if (finite['schema'], defect['schema'], primal['schema'], target['schema'],
            generators['schema'], defect['metadata'], primal['metadata'],
            defect['finite_action_sha256'], primal['finite_action_sha256'],
            primal['finite_defect_sha256'], primal['pilot_sha256'],
            metadata['target_sha256'], metadata['generators_sha256'],
            metadata['input_sha256']['inverse-residual-pilot.json.zlib'],
            generators['metadata']['target_packet_sha256'],
            generators['support_action_sha256'], metadata['N'], metadata['M'],
            metadata['last_mode'], finite['dual_first'], finite['dual_last']) != (
            'ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1', 'ACTUAL_C27_FINITE_DUAL_DEFECT_V1',
            'ACTUAL_C27_FULL_GAMMA_PRIMAL_TAIL_V1', 'ACTUAL_C27_STRUCTURED_TARGET_V1',
            'ACTUAL_C27_RESIDUAL_GENERATOR_V1', metadata, metadata,
            finite_hash, finite_hash, defect_hash, pilot_hash, target_hash,
            generator_hash, pilot_hash, target_hash, metadata['support_sha256'],
            N, M, 2048, 2049, M):
        raise ValueError('Require the same accepted common source, dual action and tail bound.')
    mw = matrix_in(defect['effective_dual_action'], defect['grid_bits'])
    residual = matrix_in(defect['residual_profile'], defect['grid_bits'])
    GO = matrix_in(generators['GO'], generators['grid_bits'])
    GP = arb_mat([[interval_ball(x, target['grid_bits']) for x in row[:2]]
                  for row in target['generators']])
    if (mw.nrows(), mw.ncols(), residual.nrows(), residual.ncols(),
            GO.nrows(), GO.ncols(), GP.nrows(), GP.ncols()) != (
            M-2048, 1, M-2048, 1, M-2048, 2, N+1, 2):
        raise ValueError('Incomplete accepted action or residual generators.')
    w = [arb(fmpq(int(x), 1 << finite['dual_bits'])) for x in finite['dual_integers']]
    weights = [arb(fmpq(int(x), 1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    if len(w) != M-2048 or len(weights) != M-N or not all(x > 0 for x in weights):
        raise ValueError('Incomplete dual coordinates or original positive weights.')
    weighted = arb_mat(M-N, 1)
    finite_norm, fixed_direction_projection = arb(0), arb(0)
    for i, m in enumerate(range(2049, M+1)):
        weight = weights[m-N-1]
        value = weight*mw[i, 0]-w[i]
        weighted[m-N-1, 0] = value
        finite_norm += weight*mw[i, 0]**2
        fixed_direction_projection += value*residual[i, 0]
    kap = interval_ball(target['kappa'], target['grid_bits'])
    projection, geometry = cauchy_projection(
        weighted, lambda m, j: GO[m-2049, j], GP,
        list(range(2049, M+1)), list(range(N+1)), kap, N, args.order)
    column = projection.transpose()
    # Reuse the saved exact v solely as a consistency check on the NEW
    # all-coordinate projection. Its original residual is not recomputed.
    direction, direction_hash = read_json(args.direction)
    if (direction['schema'], direction_hash) != (
            'ACTUAL_C27_JOINT_GAMMA_DIRECTION_V1', metadata['direction_sha256']):
        raise ValueError('Wrong accepted retained direction for the projection check.')
    v = [arb(fmpq(int(x), 1 << direction['direction_dyadic_bits']))
         for x in direction['direction_dyadic_integers']]
    if len(v) != N+1:
        raise ValueError('Incomplete original exact retained direction.')
    check = sum((column[n, 0]*v[n] for n in range(N+1)), arb(0))
    if not check.overlaps(fixed_direction_projection):
        raise ArithmeticError('Common column disagrees with the accepted same-direction action.')
    tail_norm = interval_ball(primal['bounds']['effective_action_tail_norm_upper'], primal['grid_bits'])
    far_norm_upper = (tail_norm.upper()**2).upper()
    h = interval_ball(defect['denominator'], defect['grid_bits'])
    delta_upper = (finite_norm+far_norm_upper-h).upper()
    if not delta_upper > 0:
        raise ArithmeticError('This common rank-one completion requires a positive denominator upper.')
    fields = dict(middle_action_energy=finite_norm,
                  far_action_energy_upper=far_norm_upper,
                  rank_one_delta_upper=delta_upper,
                  saved_direction_column_check=check,
                  saved_direction_profile_check=fixed_direction_projection)
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1', metadata=metadata,
                  finite_action_sha256=finite_hash, finite_defect_sha256=defect_hash,
                  primal_tail_sha256=primal_hash, target_sha256=target_hash,
                  generators_sha256=generator_hash, pilot_sha256=pilot_hash,
                  direction_sha256=direction_hash, grid_bits=args.output_bits,
                  working_bits=args.precision,
                  column=matrix_out(column, args.output_bits), finite_projection_geometry=geometry,
                  bounds={key: envelope(value, args.output_bits) for key, value in fields.items()},
                  conditions=primal['conditions']+[
                      'The same residual generators represent the full retained-to-J map R_J.',
                      'The existing finite Cauchy projection and its absolute remainder enclose every retained coordinate.'],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started,
                  scope='Only the missing common1281-coordinate column R_J^*(W_J M_J w-w) and scalar denominator allowance. '
                    'The far cross action is bounded separately; no full retained positivity, '
                    'new generic algorithm/theorem, Lean verification, cofinal positivity or RH claim.')
    write_packet(args.out, result)
    print('common column enclosed', column.nrows(), 'delta upper', delta_upper.str(24),
          'seconds', time.monotonic()-started, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(finite='full-gamma-dual-32768-finite.json.zlib',
                             defect='full-gamma-dual-32768-defect.json.zlib',
                             primal='full-gamma-primal-tail.json',
                             target='structured-joint-2048-target.json.zlib',
                             generators='joint-residual-generators-2048.json.zlib',
                             pilot='inverse-residual-pilot.json.zlib',
                             direction='joint-gamma-6765-direction.json').items():
        parser.add_argument('--'+name, type=Path, default=HERE/default)
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--order', type=int, default=512)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.output_bits != 640 or args.precision < 832 or args.order < 512:
        parser.error('Use a fresh output, the inherited640-bit grid,192 guard bits and at least512 existing Cauchy terms.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
