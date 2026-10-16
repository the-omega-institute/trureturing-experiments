"""Apply the existing structured certificate to a common rank-one target.

Read the old complete lower target and the newly enclosed common column.
Only the signed low-rank update changes; existing displacement LDL supplies
the full coupled sign test. No old producer or inverse solve is called.
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
from joint_gamma_boundary import matrix_in, matrix_out, original_target, write_packet
from structured_inverse_residual_interval import eliminate, envelope


def consume(args):
    started = time.monotonic()
    old, old_hash, metadata_bits = original_target(args)
    common, common_hash = read_json(args.common)
    rho = fmpq(args.rho)
    if rho <= 0:
        raise ValueError('Use a positive common Young parameter.')
    metadata = common['metadata']
    if (old['schema'], common['schema'], common['target_sha256'],
            metadata['target_sha256'], metadata['N'], metadata['M'],
            metadata['last_mode'], metadata['basis']) != (
            'ACTUAL_C27_STRUCTURED_TARGET_V1', 'ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1',
            old_hash, old_hash, N, M, 2048, old['metadata']['basis']):
        raise ValueError('Require the matching original complete target and common column.')
    if metadata['input_sha256'] != old['metadata']['input_sha256']:
        raise ValueError('Different actual scientific inputs.')
    for name, digest in metadata['input_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Actual input pin mismatch: '+name)
    column = matrix_in(common['column'], common['grid_bits'])
    G = matrix_in(old['generators'], old['grid_bits'])
    K = matrix_in(old['kernel'], old['grid_bits'])
    rank = K.nrows()
    if (column.nrows(), column.ncols(), G.nrows(), G.ncols(), K.ncols()) != (
            N+1, 1, N+1, 4+2*rank, rank):
        raise ValueError('Incomplete all-coordinate common update or old target.')
    # Choose exact upward rational endpoints as a single common denominator.
    bits = common['grid_bits']
    delta_up = fmpq(int(common['bounds']['rank_one_delta_upper'][1]), 1 << bits)
    tail_up = fmpq(int(common['bounds']['far_action_energy_upper'][1]), 1 << bits)
    denominator = delta_up+tail_up/rho
    if denominator <= 0 or tail_up < 0:
        raise ValueError('Positive common denominator and nonnegative far charge required.')
    inverse_denominator = arb(1)/arb(denominator)
    allowance = interval_ball(old['metadata']['complete_far_operator_allowance'],
                              metadata_bits).upper()
    P = arb_mat([[G[n, 4+j] for j in range(rank)] for n in range(N+1)])
    PK = P*K
    diagonal = []
    for n, pair in enumerate(old['diagonal']):
        far = sum((PK[n, j]*P[n, j] for j in range(rank)), arb(0))+allowance
        diagonal.append(interval_ball(pair, old['grid_bits'])-arb(rho)*far
                        +column[n, 0]**2*inverse_denominator)
    kernel = arb_mat(rank+1, rank+1)
    for i in range(rank):
        for j in range(rank):
            kernel[i, j] = (1+arb(rho))*K[i, j]
    kernel[rank, rank] = -inverse_denominator
    features = arb_mat([[P[n, j] for j in range(rank)]+[column[n, 0]]
                        for n in range(N+1)])
    generators = arb_mat([[G[n, j] for j in range(4)]
                          +[features[n, j] for j in range(rank+1)]
                          +[n*n*features[n, j] for j in range(rank+1)]
                          for n in range(N+1)])
    target_metadata = dict(old['metadata'],
                           original_target_sha256=old_hash,
                           common_column_sha256=common_hash,
                           common_young_parameter=str(rho),
                           original_metadata_grid_bits=metadata_bits,
                           original_pivot_result_sha256=hashlib.sha256(args.pivots.read_bytes()).hexdigest(),
                           reused_far_allowance=envelope(allowance, args.output_bits),
                           reused_far_allowance_grid_bits=args.output_bits,
                           rank_one_denominator=str(denominator),
                           signed_update_rank=rank+1,
                           total_generator_columns=generators.ncols())
    packet = dict(schema='ACTUAL_C27_STRUCTURED_TARGET_V1',
                  grid_bits=args.output_bits, metadata=target_metadata,
                  producer_working_bits=args.precision,
                  producer_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  diagonal=[envelope(x, args.output_bits) for x in diagonal],
                  generators=matrix_out(generators, args.output_bits),
                  kernel=matrix_out(kernel, args.output_bits),
                  kappa=envelope(interval_ball(old['kappa'], old['grid_bits']), args.output_bits))
    target_hash = write_packet(args.target_out, packet)
    # Certify the SAME outward target that is archived for later consumers.
    diagonal = [interval_ball(x, args.output_bits) for x in packet['diagonal']]
    generators = matrix_in(packet['generators'], args.output_bits)
    kernel = matrix_in(packet['kernel'], args.output_bits)
    kap = interval_ball(packet['kappa'], args.output_bits)
    verdict, failure, pivots, order = eliminate(
        diagonal, generators, kernel, kap, [n*n for n in range(N+1)],
        True, 'max-diagonal')
    result = dict(schema='ACTUAL_C27_STRUCTURED_JOINT_INTERVAL_V1',
                  N=N, M=M, last_mode=2048, basis=metadata['basis'],
                  target_packet_sha256=target_hash,
                  original_target_packet_sha256=old_hash,
                  common_column_sha256=common_hash,
                  common_young_parameter=str(rho),
                  rank_one_denominator=str(denominator),
                  working_bits=args.precision, grid_bits=args.output_bits,
                  verdict=verdict, failure_mode=failure, pivot_count=len(pivots),
                  pivot_modes=order,
                  pivot_intervals=[envelope(x, args.output_bits) for x in pivots],
                  elapsed_seconds=time.monotonic()-started,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  conditions=common['conditions']+[
                      'L0 is the matching original coefficient-uniform lower target.',
                      'The original signed far core plus its allowance majorizes the same actual residual Gram.',
                      'L0-rho*Fup+b_fin*b_fin^*/(delta_up+E_M/rho) is the common Young lower target.'],
                  scope='A full coupled structured sign test of one common rank-one/Young lower target. '
                    'A negative or inconclusive pivot does not decide the actual Schur form. '
                    'No new generic algorithm/theorem, Lean verification, cofinal positivity or RH claim.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(verdict, 'pivot count', len(pivots), 'failure', failure,
          'seconds', time.monotonic()-started, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target', type=Path, default=HERE/'structured-joint-2048-target.json.zlib')
    parser.add_argument('--pivots', type=Path, default=HERE/'structured-joint-2048-pivots.json')
    parser.add_argument('--common', type=Path, default=HERE/'full-gamma-common-column.json.zlib')
    parser.add_argument('--rho', default='1')
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--target-out', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if (args.target_out.exists() or args.out.exists() or args.precision < 832
            or args.output_bits != 640):
        parser.error('Use fresh outputs, the inherited640-bit grid and192 guard bits.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
