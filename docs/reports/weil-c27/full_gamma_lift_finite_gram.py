"""Supply the missing common finite Gram on the saved eighteen-column lift.

Reuse the residual generators, exact implicit lift, weights, channel map
and inverse. Existing certified disjoint Cauchy actions produce only the
new lift profiles and their complete weighted Gram. No selected inverse,
source action, channel/scalar supplier or full target is regenerated.
"""
import argparse
import hashlib
import platform
from pathlib import Path
import time

import flint
from flint import arb, arb_mat, ctx, fmpq

from arb_inverse_residual_obstruction import HERE, M, N
from diagonal_floor_witness import interval_ball
from full_gamma_correlated_sign import symmetric
from full_gamma_dual_middle import finite_right_action
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope


def produce(args):
    started = time.monotonic()
    helper_names = ('full_gamma_lift_finite_gram.py', 'full_gamma_dual_middle.py',
                    'full_gamma_dual.py', 'full_gamma_correlated_sign.py',
                    'gamma_channel_resolution.py', 'joint_gamma_boundary.py',
                    'structured_inverse_residual_interval.py', 'arb_arithmetic.py')
    helper_hashes = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                     for name in helper_names}
    names = ('old', 'generators', 'pilot', 'base', 'boundary', 'actions',
             'core', 'projection', 'gap')
    packets, hashes = {}, {}
    for name in names:
        packets[name], hashes[name] = read_json(getattr(args, name))
    old, generators, pilot, base, boundary, actions, core, projection, gap = (
        packets[name] for name in names)
    schemas = dict(old='ACTUAL_C27_STRUCTURED_TARGET_V1',
                   generators='ACTUAL_C27_RESIDUAL_GENERATOR_V1',
                   pilot='ACTUAL_C27_INVERSE_RESIDUAL_PILOT_V1',
                   base='ACTUAL_C27_STRUCTURED_TARGET_V1',
                   boundary='ACTUAL_C27_JOINT_GAMMA_ACTION_V1',
                   actions='ACTUAL_C27_SCHUR_LIFT_V1',
                   core='ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1',
                   projection='ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
                   gap='ACTUAL_C27_FULL_GAMMA_BASELINE_GAP_V1')
    if any(packets[name]['schema'] != schema for name, schema in schemas.items()):
        raise ValueError('Require the accepted same-residual lift and nineteen-channel packets.')
    if (generators['metadata']['target_packet_sha256'],
            base['metadata']['original_target_sha256'],
            base['metadata']['joint_boundary_action_sha256'],
            base['metadata']['common_young_parameter'],
            boundary['target_packet_sha256'], boundary['generator_packet_sha256'],
            actions['target_packet_sha256'], core['target_packet_sha256'],
            projection['original_target_sha256'], projection['base_target_sha256'],
            projection['base_core_sha256'], projection['base_lift_sha256'],
            projection['boundary_sha256']) != (
            hashes['old'], hashes['old'], hashes['boundary'], '1/32',
            hashes['old'], hashes['generators'], hashes['base'], hashes['base'],
            hashes['old'], hashes['base'], hashes['core'], hashes['actions'], hashes['boundary']):
        raise ValueError('Different residual, exact implicit lift, channel basis or base.')
    if (actions['pivot_result_sha256'], core['pivot_result_sha256']) != (
            projection['base_pivots_sha256'], projection['base_pivots_sha256']):
        raise ValueError('The implicit lift and small core must use identical accepted pivots.')
    expected = old['metadata']['input_sha256']
    if hashes['pilot'] != expected['inverse-residual-pilot.json.zlib']:
        raise ValueError('Different original exact diagonal weights.')
    if any(packets[name]['input_sha256'] != expected for name in (
            'boundary', 'core', 'projection', 'gap')) or base['metadata']['input_sha256'] != expected:
        raise ValueError('Different actual scientific inputs.')
    for name, digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Scientific input hash mismatch: '+name)
    if any(gap['packet_sha256'][name] != hashes[name] for name in (
            'old', 'base', 'projection', 'core')):
        raise ValueError('The reused far applicability result belongs to a different comparison.')
    if gap['verdict'] != 'same-Wbar-stronger-base-fixed-dual-update-excluded':
        raise ValueError('Require the accepted far applicability result, without recomputing it.')
    if int(gap['scalar_bounds']['shifted_far_quadratic_coefficient_lower'][0]) <= 0:
        raise ValueError('The archived directed far applicability lower must be strictly positive.')
    prefix, modes = actions['positive_modes'], actions['remaining_modes']
    if (len(prefix), len(modes), modes[0], sorted(prefix+modes),
            core['remaining_modes'], projection['remaining_modes'],
            old['metadata']['last_mode'], pilot['N'], pilot['M']) != (
            1263, 18, 0, list(range(N+1)), modes, modes, 2048, N, M):
        raise ValueError('Incomplete or differently ordered eighteen-column lift or finite support.')
    if (base['metadata']['original_metadata_grid_bits'],
            projection['original_metadata_grid_bits']) != (512, 512):
        raise ValueError('Keep the original allowance grid distinct from the exported input grid.')
    if (old['metadata']['N'], old['metadata']['M'], boundary['N'], boundary['M'],
            boundary['basis'], core['basis'], projection['basis'], len(boundary['bins'])) != (
            N, M, N, M, old['metadata']['basis'], old['metadata']['basis'],
            old['metadata']['basis'], 19):
        raise ValueError('Require the same complete unphased basis and original nineteen channels.')
    lift = matrix_in(actions['inverse_action'], actions['grid_bits'])
    GP = arb_mat([[interval_ball(x, old['grid_bits']) for x in row[:2]]
                  for row in old['generators']])
    GO = matrix_in(generators['GO'], generators['grid_bits'])
    beta = matrix_in(projection['channel_map_on_lift'], projection['grid_bits'])
    if (lift.nrows(), lift.ncols(), GP.nrows(), GP.ncols(), GO.nrows(), GO.ncols(),
            beta.nrows(), beta.ncols(), boundary['bins'][0], boundary['bins'][-1]) != (
            1263, 18, N+1, 2, M-2048, 2, 19, 18, [0, 1], [4181, 6765]):
        raise ValueError('Incomplete residual, original nineteen channels or exact lift.')
    V = arb_mat(N+1, 18)
    for i, n in enumerate(prefix):
        for j in range(18):
            V[n, j] = -lift[i, j]
    for j, n in enumerate(modes):
        V[n, j] = 1
    weights = [arb(fmpq(int(x), 1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    if len(weights) != M-N or not all(value > 0 for value in weights):
        raise ValueError('Require the complete positive original finite diagonal weights.')
    outside, columns, order = list(range(2049, M+1)), list(range(N+1)), 512
    kap = interval_ball(old['kappa'], old['grid_bits'])
    if not kap > 0:
        raise ValueError('The original Cauchy coefficient must be strictly positive.')
    R = arb_mat(len(outside), 18)
    profile_widths = []
    direct_checks = []
    for j in range(18):
        coefficient = arb_mat([[V[n, j]] for n in columns])
        column = finite_right_action(outside, lambda m, k: GO[m-2049, k],
                                     GP, coefficient, columns, N, order, kap)
        for i in range(len(outside)):
            R[i, j] = column[i, 0]
        profile_widths.append(max(value.upper()-value.lower() for value in column.entries()))
        if j in (0, 17):
            for m in (2049, 16384, M):
                direct = sum((kap*(GO[m-2049, 0]*GP[n, 1]-GO[m-2049, 1]*GP[n, 0])
                              /(m*m-n*n)*V[n, j] for n in columns), arb(0))
                overlap = R[m-2049, j].overlaps(direct)
                if not overlap:
                    raise ArithmeticError('The new finite lift action disagrees with a direct generator row.')
                direct_checks.append(dict(mode=m, lift_column=j, overlap=True))
        print('new lift residual column', j, 'seconds', time.monotonic()-started, flush=True)
    weighted = arb_mat([[weights[m-N-1]*R[i, j] for j in range(18)]
                        for i, m in enumerate(outside)])
    finite_gram = symmetric(R.transpose()*weighted)
    kernel = matrix_in(base['kernel'], base['grid_bits'])
    rank = kernel.nrows()-19
    if rank != len(old['kernel']) or kernel.ncols() != kernel.nrows():
        raise ValueError('Different complete far and original nineteen-channel kernel.')
    if any(not kernel[i, rank+j].is_zero() or not kernel[rank+j, i].is_zero()
           for i in range(rank) for j in range(19)):
        raise ValueError('Require the accepted independent far and channel kernel blocks.')
    inverse = symmetric(arb_mat([[-kernel[rank+i, rank+j] for j in range(19)]
                                 for i in range(19)]))
    correction = symmetric(beta.transpose()*inverse*beta)
    SJ = symmetric(finite_gram-correction)
    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], core['grid_bits']))
    HV = symmetric(matrix_in(projection['lift_gram'], projection['grid_bits']))
    far = symmetric(matrix_in(projection['far_gram_upper'], projection['grid_bits']))
    eps = fmpq(projection['exact_far_allowance'])
    ebits = base['metadata']['reused_far_allowance_grid_bits']
    if eps != fmpq(int(base['metadata']['reused_far_allowance'][1]), 1 << ebits):
        raise ValueError('Different same-source full far allowance.')
    kappa = arb(fmpq(33, 32))
    Xi = symmetric(sigma+2*kappa*arb(eps)*HV+kappa*far)
    upper = symmetric(Xi+SJ)
    matrices = dict(finite_residual_gram=finite_gram,
                    finite_channel_correction=correction,
                    channel_corrected_finite_gram=SJ,
                    far_shifted_core=Xi,
                    threshold_upper=upper)
    if any(not value.is_finite() for matrix in matrices.values() for value in matrix.entries()):
        raise ArithmeticError('The chosen finite common-Gram evaluation is not finite.')
    widths = {name: max(value.upper()-value.lower() for value in matrix.entries())
              for name, matrix in matrices.items()}
    result = dict(schema='ACTUAL_C27_LIFT_FINITE_GRAM_V1', packet_sha256=hashes,
                  input_sha256=expected, helper_sha256=helper_hashes,
                  source_sha256=helper_hashes[Path(__file__).name],
                  N=N, M=M, finite_first=2049, finite_last=M,
                  remaining_modes=modes, grid_bits=args.output_bits,
                  working_bits=args.precision, finite_cauchy_order=order,
                  direct_checks=direct_checks,
                  runtime=dict(python=platform.python_version(), flint=flint.__version__,
                               machine=platform.machine()),
                  matrices={name: matrix_out(matrix, args.output_bits)
                            for name, matrix in matrices.items()},
                  maximum_entry_width={name: envelope(value, args.output_bits)
                                       for name, value in widths.items()},
                  finite_profile_max_width=[envelope(value, args.output_bits)
                                            for value in profile_widths],
                  new_finite_lift_columns=18, new_inverse_actions=0,
                  new_channel_or_scalar_rows=0, old_source_productions_replayed=0,
                  full_target_rebuilt=False, numerical_threshold_sign='not tested',
                  verdict='same-lift-common-finite-gram-enclosed',
                  elapsed_seconds=time.monotonic()-started,
                  conditions=[
                      'Accepted residual generators describe the same exact selected-inverse residual, with exact cancellation on E.',
                      'The saved implicit lift encloses the exact same V, not a rounded replacement; all common products refer to this V.',
                      'Original positive finite W, original19channels and Beta=U*WR use the same support, normalization and weights.',
                      'The cached channel inverse is the same A_rho^-1 used by the accepted base.',
                      'Accepted complete far Gram and padding comparison are unchanged; accepted Amax>0 supplies applicability without a new sign test.',
                      'Finite M^-1<=Wc and Theta=R*WcR-R*M^-1R imply 0<=Theta<=SJ. Thus Xi+SJ is a Loewner upper, not lower, for Xi+Theta.',
                      'Matrix entries are interval enclosures of the exact common matrices, not entrywise Loewner endpoint bounds.',
                      'Interpreting a threshold result for actual whole-J sources still needs the inherited full common-domain, inverse/comparison and finite far-energy premises.'],
                  scope='New thin common finite Gram only, using existing Cauchy actions and accepted data. '
                    'No large/block inverse, new source, channel/scalar supplier, old source/action replay, '
                    'full target reproduction or threshold sign. A positive upper matrix does not prove '
                    'positivity of the actual threshold. Widths describe this evaluation, not intrinsic '
                    'conditioning. No generic novelty, Lean, actual retained/cofinal positivity, Robin-tail or RH result.')
    for name in names:
        if hashlib.sha256(getattr(args, name).read_bytes()).hexdigest() != hashes[name]:
            raise ValueError('A consumed accepted packet changed: '+name)
    for name, digest in helper_hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('A consumed helper changed: '+name)
    write_packet(args.out, result)
    for name, width in widths.items():
        print(name, 'max entry width', width.str(16), flush=True)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    defaults = dict(old='structured-joint-2048-target.json.zlib',
                    generators='joint-residual-generators-2048.json.zlib',
                    pilot='inverse-residual-pilot.json.zlib',
                    base='joint-gamma-6765-rho32-v2-target.json.zlib',
                    boundary='joint-gamma-6765-finite.json.zlib',
                    actions='joint-gamma-6765-rho32-v2-actions.json.zlib',
                    core='joint-gamma-6765-rho32-v2-schur.json',
                    projection='full-gamma-lift-projection.json.zlib',
                    gap='full-gamma-baseline-gap.json.zlib')
    for name, default in defaults.items():
        parser.add_argument('--'+name.replace('_', '-'), type=Path, default=HERE/default)
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib' or args.output_bits != 640 or args.precision < 832:
        parser.error('Use a fresh .zlib output, the inherited640-bit grid and at least192 guardbits.')
    ctx.prec = args.precision
    produce(args)


if __name__ == '__main__':
    main()
