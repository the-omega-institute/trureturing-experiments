"""Sharpen the cached threshold upper with two archived exact sources.

Reuse both source actions and the complete finite Gram. Supply only the
missing common adjoint projection and one unexported self contraction,
then use the existing inverse variational inequality and small Arb solve.
No source selection, forward lift profiles or large inverse is computed.
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
from gamma_boundary_update import read_json
from gamma_channel_resolution import cauchy_projection
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope


def consume(args):
    started = time.monotonic()
    gram, gram_hash = read_json(args.gram)
    if gram_hash != '0d28fdab5adb68146ae44735809cbfeb2cd0de55bac91143bd4548d5d69501c7':
        raise ValueError('Require the accepted exact-lift finite-Gram packet.')
    names = dict(old='structured-joint-2048-target.json.zlib',
                 generators='joint-residual-generators-2048.json.zlib',
                 pilot='inverse-residual-pilot.json.zlib',
                 base='joint-gamma-6765-rho32-v2-target.json.zlib',
                 boundary='joint-gamma-6765-finite.json.zlib',
                 actions='joint-gamma-6765-rho32-v2-actions.json.zlib',
                 core='joint-gamma-6765-rho32-v2-schur.json',
                 projection='full-gamma-lift-projection.json.zlib',
                 gap='full-gamma-baseline-gap.json.zlib')
    if set(names) != set(gram['packet_sha256']):
        raise ValueError('Incomplete accepted finite-Gram dependency map.')
    guards = {args.gram: gram_hash}
    hashes = dict(gram['packet_sha256'])
    for name, filename in names.items():
        guards[HERE/filename] = hashes[name]
    extra = dict(finite=('full-gamma-dual-32768-finite.json.zlib',
                         '4085fa7048bafc987c75b2340b56d98d556534ae74d304e61469d4e6e52af1b5'),
                 defect=('full-gamma-dual-32768-defect.json.zlib',
                         '6f320576338562a4e7338906ca66b25ab49e4a1d4460e87e26165cd9579a1154'),
                 sources=('full-gamma-two-source.json.zlib',
                          '34fca6954cdd706a11479cb68fa30995dd7586e0207c03ee6c3b68742d0896ed'))
    for name, (filename, digest) in extra.items():
        names[name], hashes[name], guards[HERE/filename] = filename, digest, digest
    for key in ('input_sha256', 'helper_sha256'):
        for filename, digest in gram[key].items():
            path = HERE/filename
            if path in guards and guards[path] != digest:
                raise ValueError('Inconsistent inherited dependency: '+filename)
            guards[path] = digest
    source = Path(__file__).resolve()
    for filename in ('arb_inverse_residual_obstruction.py', 'diagonal_floor_witness.py',
                     'full_gamma_correlated_sign.py', 'gamma_boundary_update.py',
                     'gamma_channel_resolution.py', 'joint_gamma_boundary.py',
                     'structured_inverse_residual_interval.py'):
        path = HERE/filename
        guards.setdefault(path, hashlib.sha256(path.read_bytes()).hexdigest())
    guards[source] = hashlib.sha256(source.read_bytes()).hexdigest()
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Accepted input/helper hash mismatch: '+path.name)
    # Reading accepted packets does not invoke any of their producers.
    data = {name: read_json(HERE/names[name])[0] for name in (
        'old', 'generators', 'base', 'boundary', 'actions', 'core',
        'projection', 'gap', 'finite', 'defect', 'sources')}
    old, generators, base, boundary, actions, core, projection, gap, finite, defect, sources = (
        data[name] for name in data)
    if (gram['schema'], gram['verdict'], gram['N'], gram['M'], gram['finite_first'],
            gram['finite_last'], gram['grid_bits'], finite['schema'], defect['schema'],
            sources['schema'], finite['dual_bits'], sources['new_dual_bits']) != (
            'ACTUAL_C27_LIFT_FINITE_GRAM_V1', 'same-lift-common-finite-gram-enclosed',
            N, M, 2049, M, 640, 'ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1',
            'ACTUAL_C27_FINITE_DUAL_DEFECT_V1', 'ACTUAL_C27_TWO_SOURCE_JOINT_LOWER_V1',
            384, 512):
        raise ValueError('Different accepted support, sources, entry grid or finite Gram.')
    if N != 1280 or M != 32768:
        raise ValueError('Keep the accepted finite support fixed.')
    metadata = finite['metadata']
    if (defect['metadata'] != metadata or sources['metadata'] != metadata or
            metadata['input_sha256'] != gram['input_sha256'] or
            sources['input_sha256'] != gram['input_sha256'] or
            metadata['target_sha256'] != hashes['old'] or
            metadata['generators_sha256'] != hashes['generators'] or
            metadata['basis'] != old['metadata']['basis'] or
            (metadata['N'], metadata['M'], metadata['last_mode']) != (N, M, 2048)):
        raise ValueError('Different actual source, residual, basis or scientific inputs.')
    if any((packet['dual_first'], packet['dual_last'], packet['grid_bits']) != (2049, M, 640)
           for packet in (finite, defect, sources)):
        raise ValueError('Different exact source/action support or export grid.')
    if defect['finite_action_sha256'] != hashes['finite'] or any(
            sources['packet_sha256'][name] != hashes[name] for name in (
                'finite', 'defect', 'base', 'boundary', 'actions', 'core', 'projection', 'gap')):
        raise ValueError('Different archived source actions or common lift.')
    prefix, modes = actions['positive_modes'], actions['remaining_modes']
    if (len(prefix), len(modes), modes[0], sorted(prefix+modes), gram['remaining_modes'],
            core['remaining_modes'], projection['remaining_modes'], sources['remaining_modes'],
            boundary['basis'], core['basis'], projection['basis'], len(boundary['bins']),
            base['metadata']['common_young_parameter'], gap['verdict']) != (
            1263, 18, 0, list(range(N+1)), modes, modes, modes, modes,
            metadata['basis'], metadata['basis'], metadata['basis'], 19, '1/32',
            'same-Wbar-stronger-base-fixed-dual-update-excluded'):
        raise ValueError('Different lift order, channel basis or far applicability.')
    if (projection['base_lift_sha256'], projection['base_core_sha256'],
            projection['base_target_sha256'], projection['boundary_sha256'],
            actions['target_packet_sha256'], core['target_packet_sha256']) != (
            hashes['actions'], hashes['core'], hashes['base'], hashes['boundary'],
            hashes['base'], hashes['base']):
        raise ValueError('The implicit lift and finite Gram must describe the same object.')
    if (projection['original_metadata_grid_bits'],
            base['metadata']['original_metadata_grid_bits']) != (512, 512):
        raise ValueError('Keep the allowance metadata and entry grids distinct.')

    integers = [list(map(int, finite['dual_integers'])),
                list(map(int, sources['new_dual_integers']))]
    if any(len(column) != M-2048 for column in integers):
        raise ValueError('Require both complete exact archived sources.')
    source_hashes = [hashlib.sha256(('\n'.join(map(str, column))+'\n').encode()).hexdigest()
                     for column in integers]
    if source_hashes[1] != sources['new_source_integer_sha256']:
        raise ValueError('Different archived second exact source.')
    first, second = sources['independent_source_indices']
    if integers[0][first]*integers[1][second] == integers[0][second]*integers[1][first]:
        raise ValueError('The archived exact sources must be independent.')
    source_bits = [finite['dual_bits'], sources['new_dual_bits']]
    Z = [[arb(fmpq(x, 1 << bits)) for x in column]
         for column, bits in zip(integers, source_bits)]
    action = matrix_in(sources['new_effective_action'], sources['grid_bits'])
    if (action.nrows(), action.ncols()) != (M-2048, 1):
        raise ValueError('Require the complete archived second effective action.')
    a00 = interval_ball(defect['denominator'], defect['grid_bits'])
    cross01 = interval_ball(sources['scalar_bounds']['selfadjoint_cross_forward'], 640)
    cross10 = interval_ball(sources['scalar_bounds']['selfadjoint_cross_reverse'], 640)
    if not cross01.overlaps(cross10):
        raise ArithmeticError('Archived cross terms violate common selfadjointness.')
    cross = (cross01+cross10)/2
    # This self contraction was not exported by the two-source supplier.
    a11 = sum((Z[1][i]*action[i, 0] for i in range(len(Z[1]))), arb(0))
    A = arb_mat([[a00, cross], [cross, a11]])
    second_pivot = a11-cross*cross/a00 if a00 > 0 else arb.nan()
    if not a00 > 0 or not second_pivot > 0:
        raise ArithmeticError('Inconclusive two-source positive pivots; stop without a solve retry.')

    lift = matrix_in(actions['inverse_action'], actions['grid_bits'])
    GP = arb_mat([[interval_ball(x, old['grid_bits']) for x in row[:2]]
                  for row in old['generators']])
    GO = matrix_in(generators['GO'], generators['grid_bits'])
    if (lift.nrows(), lift.ncols(), GP.nrows(), GP.ncols(), GO.nrows(), GO.ncols()) != (
            1263, 18, N+1, 2, M-2048, 2):
        raise ValueError('Incomplete common residual generators or lift.')
    V = arb_mat(N+1, 18)
    for i, n in enumerate(prefix):
        for j in range(18):
            V[n, j] = -lift[i, j]
    for j, n in enumerate(modes):
        V[n, j] = 1
    outside, columns, order = list(range(2049, M+1)), list(range(N+1)), 512
    # The helper indexes m-N-1. Its argument is UNWEIGHTED Z here,
    # with the preceding E rows zero, not WZ or WMZ-Z.
    unweighted = arb_mat(M-N, 2)
    for i, m in enumerate(outside):
        for k in range(2):
            unweighted[m-N-1, k] = Z[k][i]
    kap = interval_ball(old['kappa'], old['grid_bits'])
    if not kap > 0:
        raise ValueError('Require the original positive residual Cauchy coefficient.')
    raw, geometry = cauchy_projection(unweighted, lambda m, j: GO[m-2049, j], GP,
                                      outside, columns, kap, N, order)
    direct_checks = []
    for k in range(2):
        for n in (0, 39, N):
            direct = sum((Z[k][i]*kap*(GO[i, 0]*GP[n, 1]-GO[i, 1]*GP[n, 0])
                          /(m*m-n*n) for i, m in enumerate(outside)), arb(0))
            if not raw[k, n].overlaps(direct):
                raise ArithmeticError('Common adjoint projection disagrees with a direct generator column.')
            direct_checks.append(dict(source=k, retained_mode=n, overlap=True,
                                      projected=envelope(raw[k, n], 640),
                                      direct=envelope(direct, 640)))
    B = raw*V
    print('new two-source common projection enclosed; seconds', time.monotonic()-started, flush=True)
    solved = A.solve(B, algorithm='precond')
    K = symmetric(B.transpose()*solved)
    cached = gram['matrices']['threshold_upper']
    if len(cached) != 18 or any(len(row) != 18 for row in cached) or any(
            cached[i][j] != cached[j][i] for i in range(18) for j in range(18)):
        raise ValueError('Require the complete symmetric cached upper.')
    U2 = symmetric(matrix_in(cached, 640)-K)
    matrices = dict(source_gram=A, source_on_lift=B, inverse_gram_solution=solved,
                    inverse_projection_lower=K, threshold_upper=U2)
    if any(not x.is_finite() for matrix in matrices.values() for x in matrix.entries()):
        raise ArithmeticError('Nonfinite common projection/solve evaluation; stop.')
    determinant = U2.det()
    determinant_bits = 18*640+64
    pair = envelope(determinant, determinant_bits)
    sign = (1 if pair is not None and int(pair[0]) > 0 else
            -1 if pair is not None and int(pair[1]) < 0 else 0)
    verdict = {1: 'projection-upper-positive-actual-threshold-undecided',
               -1: 'projection-upper-negative-no-strict-negative-lower-source',
               0: 'inconclusive-projection-upper-determinant'}[sign]
    result = dict(schema='ACTUAL_C27_COMMON_PROJECTION_UPPER_V1',
                  finite_gram_sha256=gram_hash, packet_sha256=hashes,
                  input_sha256=gram['input_sha256'], remaining_modes=modes,
                  source_sha256=guards[source],
                  consumed_sha256={path.name: digest for path, digest in guards.items()},
                  exact_source_integer_sha256=source_hashes, exact_source_bits=source_bits,
                  independent_source_indices=[first, second], N=N, M=M,
                  finite_first=2049, finite_last=M, grid_bits=640, working_bits=args.precision,
                  determinant_grid_bits=determinant_bits, determinant_interval=pair,
                  determinant_sign=sign, native_determinant_calls=1,
                  common_projection_calls=1, common_projection_columns=2,
                  common_projection_geometry=geometry, direct_checks=direct_checks,
                  source_gram_pivots=[envelope(a00, 640), envelope(second_pivot, 640)],
                  source_gram_reuse='A00 and both directed cross enclosures reused; only A11 contracted.',
                  matrices={name: matrix_out(matrix, 640) for name, matrix in matrices.items()},
                  maximum_entry_width={name: envelope(max(x.upper()-x.lower()
                                         for x in matrix.entries()), 640)
                                       for name, matrix in matrices.items()},
                  new_source_actions=0, new_forward_lift_profiles=0,
                  old_scientific_productions_replayed=0, full_target_rebuilt=False,
                  runtime=dict(python=platform.python_version(), flint=flint.__version__,
                               machine=platform.machine()),
                  verdict=verdict, elapsed_seconds=time.monotonic()-started,
                  conditions=list(dict.fromkeys(gram['conditions']+gap['conditions']))+[
                      'Both exact archived sources and M actions belong to the same positive finite M_J and J support.',
                      'A=Z*M_J Z, B=Z*Y and K=B*A^-1 B refer to this same Z and exact implicit eighteen-column Y=R_J V.',
                      'The inverse variational inequality gives 0<=K<=Y*M_J^-1 Y; the inherited channel comparison gives M_J^-1<=Wc.',
                      'Hence Xi<=H=Xi+Theta<=U2=Xi+SJ-K, and the accepted seventeen-positive Xi gives at least seventeen positive directions of U2.',
                      'The accepted finite source identity n_-(Csharp)=n_+(H)-17 applies unchanged.',
                      'All entry intervals enclose common exact matrix entries; they are not entry-endpoint Loewner bounds.'],
                  scope='One missing common two-source adjoint projection and unexported A11 contraction, '
                    'then a2x2 solve and one18x18determinant using existing algorithms. A positive upper '
                    'does not decide H. A negative upper direction only excludes strict negative '
                    'directions of the same finite-J/fixed18lift lower source comparison, not actual '
                    'sufficient updates generally or semidefinite endpoint cases. An inconclusive '
                    'evaluation gives no sign. No old source/action/Gram/supplier replay, source search, '
                    'forward profiles, large inverse, generic novelty, Lean, actual retained/cofinal '
                    'positivity, Robin-tail or RH result.')
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('A consumed accepted input/helper changed: '+path.name)
    write_packet(args.out, result)
    print('determinant', determinant.str(18), flush=True)
    print(verdict, 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gram', type=Path, default=HERE/'full-gamma-lift-finite-gram.json.zlib')
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib' or args.precision < 832:
        parser.error('Use a fresh .zlib output and at least192 guardbits above the fixed640-bit grid.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
