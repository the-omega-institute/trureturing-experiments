"""Test one new coupled dual space using existing certified actions.

Reuse the old exact dual and all scalar, lift, channel and far packets.
Select one exact dyadic dual from the residual of the old negative lift
coordinate. Supply its genuinely new mixed inverse actions, effective
action and common row, then apply the existing variational and inertia
criteria jointly. No new general theorem or algorithm is supplied.
"""
import argparse
import hashlib
import math
from pathlib import Path
import time

from flint import acb, arb, arb_mat, ctx, fmpq, fmpq_mat

from arb_inverse_residual_obstruction import HERE, M, N
from arb_inverse_residual_pilot import trial_integer
from diagonal_floor_witness import interval_ball
from full_gamma_correlated_sign import symmetric
from full_gamma_directional_metric import balls, midpoint, rationals
from full_gamma_dual import boundary_entry, source_inputs
from full_gamma_dual_defect import SquaredNodeAction
from full_gamma_dual_middle import finite_right_action
from gamma_boundary_update import read_json
from gamma_channel_resolution import cauchy_projection
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope


def produce(args):
    started = time.monotonic()
    helper_names = ('full_gamma_two_source.py', 'full_gamma_dual.py',
                    'full_gamma_dual_defect.py', 'full_gamma_dual_middle.py',
                    'full_gamma_correlated_sign.py', 'full_gamma_directional_metric.py',
                    'gamma_channel_resolution.py', 'joint_gamma_boundary.py',
                    'structured_inverse_residual_interval.py', 'arb_arithmetic.py')
    helper_hashes = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                     for name in helper_names}
    (arithmetic, data, pilot, old, generators, selected, _, _, GP,
     F, _, provenance) = source_inputs(args)
    names = ('finite', 'defect', 'projection', 'actions', 'core', 'base',
             'boundary', 'actual', 'certificate', 'gap')
    packets, hashes = {}, {}
    for name in names:
        packets[name], hashes[name] = read_json(getattr(args, name))
    finite, defect, projection, actions, core, base, boundary, actual, certificate, gap = (
        packets[name] for name in names)
    schemas = dict(finite='ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1',
                   defect='ACTUAL_C27_FINITE_DUAL_DEFECT_V1',
                   projection='ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
                   actions='ACTUAL_C27_SCHUR_LIFT_V1',
                   core='ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1',
                   base='ACTUAL_C27_STRUCTURED_TARGET_V1',
                   boundary='ACTUAL_C27_JOINT_GAMMA_ACTION_V1',
                   actual='ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
                   certificate='ACTUAL_C27_FULL_GAMMA_DIRECTIONAL_METRIC_V1',
                   gap='ACTUAL_C27_FULL_GAMMA_BASELINE_GAP_V1')
    if any(packets[name]['schema'] != schema for name, schema in schemas.items()):
        raise ValueError('Require accepted same-source dual, lift and comparison packets.')
    if (finite['metadata'], defect['metadata'], defect['finite_action_sha256'],
            projection['base_lift_sha256'], projection['base_core_sha256'],
            projection['base_target_sha256'], projection['boundary_sha256'],
            actual['packet_sha256']['defect'], actual['packet_sha256']['projection'],
            certificate['packet_sha256']['core'], base['metadata']['common_young_parameter'],
            actions['target_packet_sha256'], core['target_packet_sha256']) != (
            provenance, provenance, hashes['finite'], hashes['actions'], hashes['core'],
            hashes['base'], hashes['boundary'], hashes['defect'], hashes['projection'],
            hashes['core'], '1/32', hashes['base'], hashes['base']):
        raise ValueError('Different actual source, exact lift, reference or nineteen-channel base.')
    if gap['packet_sha256']['projection'] != hashes['projection'] or gap['packet_sha256']['core'] != hashes['core']:
        raise ValueError('Different accepted stronger-base padding comparison.')
    expected = old['metadata']['input_sha256']
    for packet in (projection, core, boundary, actual, certificate):
        if packet['input_sha256'] != expected:
            raise ValueError('Different scientific inputs.')
    for name, digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Scientific input hash mismatch: '+name)
    prefix, modes = actions['positive_modes'], actions['remaining_modes']
    if (len(prefix), len(modes), modes[0], sorted(prefix+modes), finite['dual_first'],
            finite['dual_last'], defect['dual_first'], defect['dual_last']) != (
            1263, 18, 0, list(range(N+1)), 2049, M, 2049, M):
        raise ValueError('Incomplete or differently ordered retained lift and dual support.')
    lift = matrix_in(actions['inverse_action'], actions['grid_bits'])
    V = arb_mat(N+1, 18)
    for i, n in enumerate(prefix):
        for j in range(18):
            V[n, j] = -lift[i, j]
    for j, n in enumerate(modes):
        V[n, j] = 1
    weights = [arb(fmpq(int(x), 1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    outside = list(range(2049, M+1))
    GO = matrix_in(generators['GO'], generators['grid_bits'])
    kap = interval_ball(old['kappa'], old['grid_bits'])
    selected_source_hash = None
    if args.reuse_source:
        selected_source, selected_source_hash = read_json(args.reuse_source)
        if (selected_source['schema'], selected_source['metadata'],
                selected_source['dual_first'], selected_source['dual_last'],
                selected_source['new_dual_bits']) != (
                'ACTUAL_C27_TWO_SOURCE_JOINT_LOWER_V1', provenance, 2049, M, 512):
            raise ValueError('Different exact source supplied for independent consumption.')
        wi = [int(x) for x in selected_source['new_dual_integers']]
        shift = selected_source['proposal_power_two_shift']
    else:
        r1 = finite_right_action(outside, lambda m, j: GO[m-2049, j], GP,
                                 arb_mat([[V[n, 0]] for n in range(N+1)]),
                                 list(range(N+1)), N, 256, kap)
        proposal = [weights[m-N-1]*r1[i, 0] for i, m in enumerate(outside)]
        amplitude = max(abs(float(x.mid())) for x in proposal)
        if not math.isfinite(amplitude) or amplitude == 0:
            raise ArithmeticError('Require a nonzero finite new-source proposal.')
        shift = -math.floor(math.log2(amplitude))
        wi = [trial_integer(x*arb(2)**shift, 512) for x in proposal]
    if len(wi) != len(outside) or not any(wi):
        raise ValueError('Require a complete nonzero new exact source.')
    w1 = [arb(fmpq(x, 1 << 512)) for x in wi]
    w0q = [fmpq(int(x), 1 << finite['dual_bits']) for x in finite['dual_integers']]
    w0 = [arb(x) for x in w0q]
    first = next(i for i, x in enumerate(w0q) if x)
    independent = next((j for j in range(len(wi))
                        if w0q[first]*wi[j]-w0q[j]*wi[first]), None)
    if independent is None:
        raise ArithmeticError('The proposed source adds no exact independent direction.')
    source_hash = hashlib.sha256(('\n'.join(map(str, wi))+'\n').encode()).hexdigest()
    print('new exact independent source', source_hash, 'seconds', time.monotonic()-started, flush=True)

    bits = finite['grid_bits']
    gp = {row['n']: (interval_ball(row['s'], bits), interval_ball(row['a'], bits))
          for row in finite['boundary_scalars']}
    periodic = {row['n']: interval_ball(row['periodic'], bits) for row in finite['boundary_scalars']}
    all_modes = list(range(N+1, M+1))
    if set(gp) != set(all_modes) or len(gp) != len(finite['boundary_scalars']):
        raise ValueError('Require the complete unchanged Gamma/pole scalar array.')
    B0 = arb_mat([[boundary_entry(arithmetic, m, n, gp, periodic) for n in selected]
                  for m in selected])
    for i, m in enumerate(selected):
        B0[i, i] += 1/weights[m-N-1]
    D = F-B0
    action = SquaredNodeAction(N+1, M)
    gamma = arithmetic.d/arithmetic.pi**2
    omega = [2*arithmetic.pi*n/arithmetic.d for n in all_modes]
    phi_gp = [-omega[i]*gp[n][0] for i, n in enumerate(all_modes)]
    diag_gp = [gp[n][1]-gp[n][0]/(arithmetic.pi*n)-periodic[n] for n in all_modes]
    phi_h = [-omega[i]*data[n][0] for i, n in enumerate(all_modes)]
    diag_h = [data[n][1]-data[n][0]/(arithmetic.pi*n) for n in all_modes]
    source = [arb(0)]*768+w1
    hw = action.apply(source, phi_h, diag_h, gamma)
    gw = action.apply(source, phi_gp, diag_gp, gamma)
    direct_modes = [selected[0], selected[-1], M]
    for m in direct_modes:
        direct_h = sum((arithmetic.even_entry(m, n, data)*w1[i]
                        for i, n in enumerate(outside)), arb(0))
        direct_g = sum((boundary_entry(arithmetic, m, n, gp, periodic)*w1[i]
                        for i, n in enumerate(outside)), arb(0))
        if not hw[m-N-1].overlaps(direct_h) or not gw[m-N-1].overlaps(direct_g):
            raise ArithmeticError('New convolution disagrees with the same actual direct row.')
    z = arb_mat([[hw[i]] for i in range(768)])
    mixed_b = arb_mat([[gw[i]] for i in range(768)])
    p, q = F.solve(z, algorithm='precond'), D.solve(z-mixed_b, algorithm='precond')
    g_action = action.apply([-x for x in q.entries()]+w1, phi_gp, diag_gp, gamma)
    h_action = action.apply(list((q-p).entries())+[arb(0)]*len(w1), phi_h, diag_h, gamma)
    u1 = [w1[i]/weights[m-N-1]+g_action[i+768]+h_action[i+768]
          for i, m in enumerate(outside)]
    u0 = list(matrix_in(defect['effective_dual_action'], defect['grid_bits']).entries())
    cross01 = sum((w0[i]*u1[i] for i in range(len(w1))), arb(0))
    cross10 = sum((w1[i]*u0[i] for i in range(len(w1))), arb(0))
    if not cross01.overlaps(cross10):
        raise ArithmeticError('New and cached mixed actions violate the common selfadjoint cross.')
    h00 = interval_ball(defect['denominator'], defect['grid_bits'])
    h11 = sum((w1[i]*u1[i] for i in range(len(w1))), arb(0))
    h = arb_mat([[h00, (cross01+cross10)/2], [(cross01+cross10)/2, h11]])
    U = arb_mat([[u0[i], u1[i]] for i in range(len(w1))])
    weighted_U = arb_mat([[weights[m-N-1]*u0[i], weights[m-N-1]*u1[i]]
                          for i, m in enumerate(outside)])
    middle_gram = symmetric(U.transpose()*weighted_U)
    weighted = arb_mat(M-N, 1)
    for i, m in enumerate(outside):
        weighted[m-N-1, 0] = weights[m-N-1]*u1[i]-w1[i]
    common_row, geometry = cauchy_projection(weighted, lambda m, j: GO[m-2049, j],
                                             GP, outside, list(range(N+1)), kap, N, 512)
    commonV = common_row*V
    print('new effective action and common row enclosed; seconds', time.monotonic()-started, flush=True)

    bins = boundary['bins']
    if len(bins) != 19 or bins[0] != [0, 1] or bins[-1][1] != 6765:
        raise ValueError('Require unchanged original nineteen channel bins.')
    factors = [((1-(-arb(fmpq(4*a+1, 2))*arithmetic.d).exp())/((b-a)*arithmetic.d)).sqrt()
               for a, b in bins]
    endpoints = sorted({value for pair in bins for value in pair})
    phi1 = arb_mat(1, 19)
    for i, m in enumerate(outside):
        im = arithmetic.pi*m/arithmetic.d
        values = {e: acb(arb(e)+arb(fmpq(1, 4)), im).digamma().real for e in endpoints}
        weight = (-1)**m*weights[m-N-1]*u1[i]
        for j, (a, b) in enumerate(bins):
            phi1[0, j] += weight*factors[j]*(values[b]-values[a])
    phi0 = matrix_in(actual['actual_channel_row'], actual['grid_bits'])
    phi = arb_mat([[phi0[0, j] for j in range(19)], [phi1[0, j] for j in range(19)]])
    kernel = matrix_in(base['kernel'], base['grid_bits'])
    rank = kernel.nrows()-19
    inv = symmetric(arb_mat([[-kernel[rank+i, rank+j] for j in range(19)] for i in range(19)]))
    betaV = matrix_in(projection['channel_map_on_lift'], projection['grid_bits'])
    old_commonV = matrix_in(projection['common_column_on_lift'], projection['grid_bits'])
    common = arb_mat([[old_commonV[0, j] for j in range(18)], [commonV[0, j] for j in range(18)]])
    B = common-phi*inv*betaV
    KJ = symmetric(middle_gram-h-phi*inv*phi.transpose())

    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], core['grid_bits']))
    lift_gram = symmetric(matrix_in(projection['lift_gram'], projection['grid_bits']))
    G = symmetric(matrix_in(projection['far_gram_upper'], projection['grid_bits']))
    epsilon = fmpq(projection['exact_far_allowance'])
    ebits = base['metadata']['reused_far_allowance_grid_bits']
    if epsilon != fmpq(int(base['metadata']['reused_far_allowance'][1]), 1 << ebits):
        raise ValueError('Different transported complete far allowance.')
    kappa = arb(fmpq(33, 32))
    sigma_max = symmetric(sigma+2*kappa*arb(epsilon)*lift_gram)
    a = sigma_max[0, 0]
    if not a < 0:
        raise ArithmeticError('The stronger compression upper needs its negative first pivot.')
    T = arb_mat([[-sigma_max[0, j+1]/a if i == 0 else arb(int(i == j+1))
                  for j in range(17)] for i in range(18)])
    P = symmetric(T.transpose()*sigma_max*T)
    H = B*T
    Rref = fmpq_mat([[fmpq(x) for x in row] for row in certificate['reference_inverse']])
    rhs = fmpq_mat([[midpoint(envelope(H[j, i], args.output_bits), args.output_bits)
                     for j in range(2)] for i in range(17)])
    Xq = Rref*rhs
    X = balls(Xq)
    B0col = arb_mat([[B[i, 0]] for i in range(2)])
    qtrial = symmetric(B0col*B0col.transpose()/a+H*X+X.transpose()*H.transpose()-X.transpose()*P*X)
    g = T*X
    for j in range(2):
        g[0, j] += B[j, 0]/a
    A = kappa-kappa**2*G[0, 0]/(-a)
    if not A > 0:
        raise ArithmeticError('Require positive same-source far quadratic coefficient.')
    penalty = symmetric(kappa**2/A*(g.transpose()*G*g))
    C = symmetric(KJ+qtrial-penalty)
    first_pivot = C[0, 0]
    second_pivot = C[1, 1]-C[1, 0]*C[0, 1]/first_pivot if not first_pivot.contains(0) else arb.nan()
    positive = first_pivot > 0 and second_pivot > 0
    print('joint lower matrix', C, flush=True)
    scalar = dict(independence_minor=arb(w0q[first]*fmpq(wi[independent], 1 << 512)
                  -w0q[independent]*fmpq(wi[first], 1 << 512)),
                  selfadjoint_cross_forward=cross01, selfadjoint_cross_reverse=cross10,
                  far_coefficient=A, first_lower_pivot=first_pivot,
                  second_lower_pivot=second_pivot)
    result = dict(schema='ACTUAL_C27_TWO_SOURCE_JOINT_LOWER_V1', packet_sha256=hashes,
                  reused_exact_source_sha256=selected_source_hash,
                  metadata=provenance, input_sha256=expected, remaining_modes=modes,
                  helper_sha256=helper_hashes, new_source_integer_sha256=source_hash,
                  grid_bits=args.output_bits, working_bits=args.precision,
                  dual_first=2049, dual_last=M, new_dual_bits=512,
                  new_dual_integers=list(map(str, wi)), proposal_power_two_shift=shift,
                  independent_source_indices=[first, independent],
                  direct_sample_modes=direct_modes,
                  new_mixed_actual=matrix_out(z, args.output_bits),
                  new_mixed_comparison=matrix_out(mixed_b, args.output_bits),
                  new_actual_inverse_action=matrix_out(p, args.output_bits),
                  new_comparison_inverse_action=matrix_out(q, args.output_bits),
                  new_effective_action=matrix_out(arb_mat([[x] for x in u1]), args.output_bits),
                  new_common_row=matrix_out(common_row, args.output_bits),
                  new_channel_row=matrix_out(phi1, args.output_bits),
                  finite_row_on_lift=matrix_out(B, args.output_bits),
                  finite_charge_matrix=matrix_out(KJ, args.output_bits),
                  variational_trial=rationals(Xq),
                  joint_lower_matrix=matrix_out(C, args.output_bits),
                  scalar_bounds={name: envelope(value, args.output_bits) for name, value in scalar.items()
                                 if value.is_finite()},
                  common_projection_geometry=geometry,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started,
                  verdict=('specified-two-source-stronger-base-update-excluded' if positive
                           else 'inconclusive-two-source-joint-lower'),
                  conditions=gap['conditions']+[
                      'The new dyadic source is an exact independent finite trial, not an exact residual.',
                      'The common finite-input domain contains this entire source span and both effective actions have finite weighted far energy.',
                      'The new p and q are certified actions on the new mixed right sides of the same actual F and D.',
                      'All cross terms use this same pair of actual sources, common form, weights and channels.',
                      'The full signed far Gram plus its entire allowance majorizes Rfar V for every coefficient.',
                      'The accepted old17positive inertia and two-sided far-padding comparison apply unchanged.'],
                  scope='Conditional numerical application of existing block inertia, variational inversion, '
                    'weighted Cauchy and square completion. A positive verdict excludes simultaneous updates '
                    'from this exact two-source space for the same19channel stronger base, every positive '
                    'denominator matrix D>=Kactual. Other sources, majorants, actual S_N sign, full retained/'
                    'cofinal positivity, Robin tail and RH are not decided. No generic novelty or Lean claim.')
    for name, digest in helper_hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('A consumed helper changed during production: '+name)
    write_packet(args.out, result)
    for name, value in scalar.items():
        print(name, value.str(18), flush=True)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    defaults = dict(target='structured-joint-2048-target.json.zlib',
                    pivots='structured-joint-2048-pivots.json',
                    support='support-2048-finite-action.json.zlib',
                    generators='joint-residual-generators-2048.json.zlib',
                    direction='joint-gamma-6765-direction.json',
                    finite='full-gamma-dual-32768-finite.json.zlib',
                    defect='full-gamma-dual-32768-defect.json.zlib',
                    projection='full-gamma-lift-projection.json.zlib',
                    actions='joint-gamma-6765-rho32-v2-actions.json.zlib',
                    core='joint-gamma-6765-rho32-v2-schur.json',
                    base='joint-gamma-6765-rho32-v2-target.json.zlib',
                    boundary='joint-gamma-6765-finite.json.zlib',
                    actual='full-gamma-actual-channel-projection.json.zlib',
                    certificate='full-gamma-directional-metric.json.zlib',
                    gap='full-gamma-baseline-gap.json.zlib')
    for name, default in defaults.items():
        parser.add_argument('--'+name, type=Path, default=HERE/default)
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--reuse-source', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib' or args.output_bits != 640 or args.precision < 832:
        parser.error('Use a fresh .zlib output and the inherited640-bit grid with192 guardbits.')
    ctx.prec = args.precision
    produce(args)


if __name__ == '__main__':
    main()
