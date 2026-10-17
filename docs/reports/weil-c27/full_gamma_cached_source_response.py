"""Transport saved finite outputs to a new SAME-source signed response.

Reuse Cauchy and old-jet samples, an explicit Cauchy-kernel interpolation,
saved real tail kernels and original physical errors. No old supplier runs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, arb_mat, arb_poly, ctx, fmpq, fmpq_poly

from diagonal_floor_witness import interval_ball
from full_gamma_coherent_split_capacity import freeze_upper
from full_gamma_joint_plane_paid_floor import plane_floor
from full_gamma_negative_index_joint import freeze_lower
from full_gamma_signed_mu import require_finite, verify_guards, write_packet
from joint_gamma_boundary import matrix_in, matrix_out
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS, M, R, ORDER = 640, 32768, 56, 13
PINS = {
    'full-gamma-signed-jet-interface.json.zlib':
        '154949d8ab3e1f2875de531aa63a9c5a3c3d7d133640a5ce8b1d23c98d67565a',
    'full-gamma-signed-mu-direct-sine.json.zlib':
        '6475e5f63dba423b2f1ee7e488581f15ef8d53de83a22ea064ac2e393b124d55',
    'full-gamma-signed-mu-direct-sine.json.finite.zlib':
        '43a29c6710cae4b4b7ad2fce9433ebba4ca0bdec9dfea968d668b13d3a469bc8',
    'full-gamma-signed-mu-direct-sine.json.cauchy.zlib':
        'dbfedf8f55f74b4fa1659ff02e57559fea52e3250379c8f99ed1807ee91245f9',
    'full-gamma-signed-correction.json.zlib':
        'd527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd',
    'full-gamma-zero-channel-pairing.json.zlib':
        'c4cbbead777f652812c0aa0bd238bf08053daa52be7afbedfa1321da19847a9b',
    'full-gamma-zero-channel-jet-refine.json.zlib':
        'ee401d13ef503a44089c33e958b8f90dab09aefa843050f09713ee6a33729861',
    'full-gamma-zero-channel-cross-refine.json.zlib':
        '5cfba16d5428c1c61618a861a8167ca6828f46142d7bd7993d4d84b20de6e980',
    'full-gamma-joint-plane-paid-floor.json.zlib':
        '3bdf5f3f81daeeda99c5fa71f2a48ee64c84c09da0f20070b45b0ad781aeec22',
    'structured-joint-2048-target.json.zlib':
        'ea4a761f76b59a11b2d04faf80531901a5164f4145c2424a4e79889ae3f72e22',
    'tail-moments.json.zlib':
        'fc91b487a48f82cb276c0823f28c4a38432f12f822d6e0ed2644a67da7c22ad1',
    'full-gamma-hankel-tail-energy.json.zlib':
        'd09a5ab254d164cfd40982b31e1a155e31d6288b705d916bf248236a0d527606',
}


def recover_symbols(jet, old_features, scale):
    """Recover the SAME approximate symbol, not the actual physical symbol."""
    if scale < 2 or len(jet) != scale or len(old_features) != 112:
        raise ValueError('Require complete saved finite jet and old short features.')
    require_finite(list(jet)+list(old_features))
    aa, bb = old_features[:56], old_features[56:]
    denominator = abs(aa[0]).lower()-sum((abs(v).upper() for v in aa[1:]), arb(0))
    if not denominator > 0:
        raise ArithmeticError('Old A polynomial has no certified nonzero denominator.')
    xs = [arb(fmpq(scale, scale+1+i)) for i in range(scale)]
    squares = [x*x for x in xs]
    av, bv = (arb_poly(h).evaluate(squares, algorithm='iter') for h in (aa, bb))
    symbols = [(x*b-arb.pi()*(scale+1+i)*(-1)**(scale+1+i)*j/2)/a
               for i, (x, a, b, j) in enumerate(zip(xs, av, bv, jet))]
    require_finite(symbols)
    return symbols, freeze_lower(denominator.lower())


def interpolate_cauchy(samples, scale, source_l1, order=ORDER):
    """Explicit Lagrange interpolation of U(x)=sum z_n/(1+(n/M)x).

    Source coefficients have support 0<n<=M and sum|z_n|<=source_l1.
    All nodes are >=1/2; the returned remainder holds on [0,1/2].
    These are inherited source premises, not inferred from the samples.
    """
    require_finite(list(samples)+[arb(source_l1)])
    if scale < order or order < 2 or len(samples) != scale or not source_l1 >= 0:
        raise ValueError('Invalid fixed-band Cauchy interpolation domain.')
    modes = [scale+1+i*(scale-1)//(order-1) for i in range(order)]
    nodes = [fmpq(scale, mode) for mode in modes]
    if len(set(map(str, nodes))) != order:
        raise ValueError('Repeated interpolation nodes.')
    polynomial, product = arb_poly(), fmpq(1)
    for i, (mode, node) in enumerate(zip(modes, nodes)):
        basis, divisor = fmpq_poly([1]), fmpq(1)
        for j, other in enumerate(nodes):
            if j != i:
                basis *= fmpq_poly([-other, 1])
                divisor *= node-other
        basis /= divisor
        polynomial += arb_poly(basis)*samples[mode-scale-1]*mode
        product *= node/(1+node)
    remainder = freeze_upper(arb(source_l1)*arb(product))
    require_finite(list(polynomial.coeffs())+[arb(remainder)])
    return polynomial, modes, nodes, product, remainder


def finite_pairing(samples, symbols, features, scale, kappa, weight):
    """New linear functional on old samples; no old jet action is repeated."""
    if (scale < 2 or len(samples) != scale or len(symbols) != scale
            or len(features) != 112 or not kappa > 0 or not weight > 0):
        raise ValueError('Invalid common finite functional domain.')
    require_finite(list(samples)+list(symbols)+list(features))
    xs = [arb(fmpq(scale, scale+1+i)) for i in range(scale)]
    squares = [x*x for x in xs]
    av, bv = (arb_poly(h).evaluate(squares, algorithm='iter')
              for h in (features[:56], features[56:]))
    total = sum((c*(x*b-s*a)/(scale+1+i)
                 for i, (c, x, a, b, s) in
                 enumerate(zip(samples, xs, av, bv, symbols))), arb(0))
    total *= arb(kappa)*arb(weight)/arb.pi()
    require_finite([total])
    return total


def cached_tail(polynomial, features, smooth, kernels, phase_cap, scale, kappa, weight):
    """Reuse K_q=M^(q-1) sum_{m>=2M+1}m^-q and saved Abel phase cap."""
    if (len(features) != 112 or not phase_cap >= 0 or scale < 2
            or not kappa > 0 or not weight > 0):
        raise ValueError('Invalid cached tail domain.')
    even, odd = [arb(0)]*111, [arb(0)]*112
    for k in range(56):
        even[2*k], odd[2*k+1] = features[k], features[56+k]
    f0 = polynomial*(arb_poly(odd)-smooth*arb_poly(even))
    f1 = polynomial*arb_poly(even)
    if any(str(p+2) not in kernels for p in range(len(f0))):
        raise ValueError('Saved real kernels do not cover the new polynomial.')
    total = sum((v*kernels[str(p+2)] for p, v in enumerate(f0.coeffs())), arb(0))
    total *= arb(kappa)*arb(weight)/(arb.pi()*scale)
    ratio = arb(fmpq(scale, 2*scale+1))
    cap = sum((abs(v).upper()*ratio**p for p, v in enumerate(f1.coeffs())), arb(0))
    oscillatory = freeze_upper(2*arb(kappa)*arb(weight)*cap*arb(phase_cap)
                               /(arb.pi()*(2*scale+1)**2))
    require_finite([total, cap, arb(oscillatory)])
    return total, cap, oscillatory


def consume(out):
    started, source = time.monotonic(), Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            invoked.add(Path(filename).resolve())
    guards = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    packets = {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted cached-source input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw))
    (jet, mu, finite, cauchy, signed, zero, refine, cross, plane, target,
     seed, hankel) = [packets[name] for name in PINS]
    inherited = dict(PINS)
    for packet in (mu, signed, zero, refine, cross, plane):
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256'):
            for name, digest in packet.get(field, {}).items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting cached-source binding: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked cached-source binding: '+name)
        guards[path] = digest
    verify_guards(guards)
    expected = ['ACTUAL_C27_SIGNED_JET_INTERFACE_V1',
        'ACTUAL_C27_SAME_SOURCE_SIGNED_MU_ABEL_V1', 'ACTUAL_C27_SIGNED_MU_FINITE_BLOCK_V1',
        'ACTUAL_C27_SIGNED_MU_CAUCHY_BLOCK_V1', 'ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1',
        'ACTUAL_C27_ZERO_CHANNEL_PAIRING_BOUND_V1', 'ACTUAL_C27_ZERO_CHANNEL_JET_REFINEMENT_V1',
        'ACTUAL_C27_ZERO_CHANNEL_CROSS_REFINEMENT_V1', 'ACTUAL_C27_JOINT_PLANE_PAID_FLOOR_V1',
        'ACTUAL_C27_STRUCTURED_TARGET_V1', 'ACTUAL_WEIL_NORMALIZED_TAIL_MOMENTS_V1',
        'ACTUAL_C27_SIGNED_HANKEL_TAIL_ENERGY_V1']
    if ([p['schema'] for p in packets.values()] != expected
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any((p['source_bits'], p['source_row'], p['source_integer_sha256']) !=
                   (384, 0, mu['source_integer_sha256'])
                   for p in (jet, finite, cauchy, signed, zero, refine, cross, plane, hankel))
            or target['metadata']['input_sha256'] != jet['input_sha256']
            or zero['input_sha256'] != jet['input_sha256']
            or (seed['c'], seed['M']) != (27, M)
            or seed['analytic_grid_bits'] != mu['symbol_metadata_grid_bits']
            or jet['original_metadata_grid_bits'] != 512
            or any(p['feature_order'] != 'A0..A55,B0..B55' for p in (jet, finite, mu, signed))
            or (target['metadata']['N'], target['metadata']['M'], target['metadata']['last_mode'])
                != (1280, M, 2048)):
        raise ValueError('Different saved source, common lift, grids or physical support.')
    for p in (finite, cauchy, mu):
        if (p['finite_first_mode'], p['finite_last_mode'], p['normalization_scale']) != (M+1, 2*M, M):
            raise ValueError('Different saved finite Cauchy/jet band.')
        if (p['first_source_mode'], p['last_source_mode'],
                p['physical_residual_last_mode']) != (1281, M, 2048):
            raise ValueError('Different original Cauchy source support.')
    if (mu['analytic_first_mode'] != 2*M+1 or mu['source_power_order'] != ORDER
            or mu['finite_checkpoint'] != dict(path='full-gamma-signed-mu-direct-sine.json.finite.zlib',
                sha256=PINS['full-gamma-signed-mu-direct-sine.json.finite.zlib'])
            or mu['cauchy_checkpoint'] != finite['cauchy_checkpoint']
            or finite['cauchy_checkpoint'] != dict(path='full-gamma-signed-mu-direct-sine.json.cauchy.zlib',
                sha256=PINS['full-gamma-signed-mu-direct-sine.json.cauchy.zlib'])
            or not finite['symbol_is_approximate'] or not finite['errors_paid_in_final_packet']):
        raise ValueError('Different finite checkpoints or approximate-symbol contract.')
    rational = lambda p, key: fmpq(*map(int, p['exact_rationals'][key]))
    ball = lambda p, key: interval_ball(p['scalar_bounds'][key], BITS)
    upper = lambda p, key: fmpq(int(p['scalar_bounds'][key][1]), 1 << BITS)
    kappa, weight = rational(mu, 'kappa'), rational(mu, 'original_far_weight')
    if kappa != fmpq(33, 32) or kappa != rational(cross, 'kappa') or not weight > 0:
        raise ValueError('Different original weight or response units.')
    samples = [interval_ball(row, BITS) for row in cauchy['finite_cauchy_source_outputs']]
    old_jet = [interval_ball(row, BITS) for row in finite['finite_approximate_jet_outputs']]
    hb = list(matrix_in(jet['short_jet_coefficients'], BITS).entries())
    hd_matrix = matrix_in(zero['new_short_jet_coefficients'], BITS)
    hd = list(hd_matrix.entries())
    left = matrix_in(signed['left_midpoint_features'], BITS)
    kernel = matrix_in(target['kernel'], BITS)
    if (len(hb) != 112 or (hd_matrix.nrows(), hd_matrix.ncols()) != (112, 1)
            or (left.nrows(), left.ncols(), kernel.nrows(), kernel.ncols()) != (112, 1, 112, 112)
            or len(seed['normalized_polynomial']) != 64):
        raise ValueError('Incomplete cached common feature arrays.')
    # All actual-data conversions/actions start below the metadata-only review interface.
    symbols, denominator = recover_symbols(old_jet, hb, M)
    z_l1 = upper(mu, 'cached_source_l1_upper')
    polynomial, modes, nodes, product, remainder = interpolate_cauchy(samples, M, z_l1)
    smooth = arb_poly([interval_ball(row, BITS) for row in seed['normalized_polynomial']])
    sigma = interval_ball(seed['uniform_symbol_error'], seed['analytic_grid_bits']).upper()
    kernels = {key: interval_ball(value, BITS) for key, value in mu['real_hurwitz_kernels'].items()}
    phase_cap = upper(mu, 'inverse_phase_gap_weight_cap')
    finite_mu = finite_pairing(samples, symbols, hd, M, kappa, weight)
    tail_mu, tail_cap, osc_error = cached_tail(polynomial, hd, smooth, kernels,
                                              phase_cap, M, kappa, weight)
    eps_d = rational(refine, 'new_sharp_forward_error_upper')
    norm_d2 = rational(cross, 'new_representer_norm2_upper')
    if not norm_d2 > 0 or not eps_d >= 0 or not sigma >= 0:
        raise ValueError('Invalid saved actual representer or symbol-error cap.')
    pcap = upper(hankel, 'principal_weighted_norm_upper')
    delta = freeze_upper(arb(weight).sqrt()*arb(remainder)/(2*arb(2*M).sqrt()))
    ratio = arb(fmpq(M, M+1))
    aa = sum((abs(v).upper()*ratio**(2*k) for k, v in enumerate(hd[:R])), arb(0))
    beta = freeze_upper(2*arb(weight).sqrt()*sigma*aa/(arb.pi()*arb(M).sqrt()))
    jcap = arb(norm_d2).sqrt()+arb(eps_d)
    functional_error = freeze_upper(arb(kappa)*(jcap*arb(delta)+(arb(pcap)+arb(delta))*arb(beta)))
    actual_proxy = finite_mu+tail_mu+arb(0, arb(functional_error+osc_error))
    # New cached physical-kernel action; no source feature reducer or old cross runs.
    left_action = left.transpose()*kernel
    mixed = arb(kappa)*(left_action*hd_matrix)[0, 0]
    old_lift = rational(jet, 'lift_norm_l1_upper')
    new_lift = rational(zero, 'new_direction_lift_norm_upper')
    if not old_lift > 0 or not new_lift >= 0:
        raise ValueError('Missing homogeneous physical-error lift domains.')
    mixed_error = freeze_upper(arb(rational(signed, 'mixed_cross_error_upper'))
                               *arb(new_lift)/arb(old_lift))
    eps_sm = upper(hankel, 'smooth_remainder_weighted_norm_upper')
    if min(pcap, eps_sm, mixed_error) < 0:
        raise ValueError('Missing nonnegative saved physical/source error caps.')
    smooth_error = freeze_upper(arb(kappa)*arb(eps_sm)*arb(norm_d2).sqrt())
    forward_error = freeze_upper(arb(kappa)*arb(pcap)*arb(eps_d))
    response_error = mixed_error+smooth_error+forward_error
    actual_zeta = actual_proxy+mixed+arb(0, arb(response_error))
    previous_lo, previous_hi = (rational(cross, 'zero_pairing_'+key) for key in ('lower', 'upper'))
    zlo = max(previous_lo, freeze_lower(actual_zeta.lower()))
    zhi = min(previous_hi, freeze_upper(actual_zeta.upper()))
    if zlo > zhi:
        raise ArithmeticError('New actual source response conflicts with saved source interval.')
    scalars = dict(finite_principal_functional=finite_mu, tail_principal_functional=tail_mu,
        tail_oscillatory_coefficient_cap=tail_cap, actual_proxy_functional=actual_proxy,
        mixed_signed_kernel_cross=mixed, actual_far_zero_pairing=actual_zeta)
    exact = dict(old_A_absolute_lower=denominator, source_l1_upper=z_l1,
        interpolation_remainder_upper=remainder, interpolation_product=product,
        source_tail_norm_error_upper=delta, symbol_jet_norm_error_upper=beta,
        functional_analytic_error_upper=functional_error, oscillatory_tail_error_upper=osc_error,
        mixed_cross_error_upper=mixed_error, smooth_error_upper=smooth_error,
        proxy_forward_error_upper=forward_error, complete_response_error_upper=response_error,
        zeta_lower=zlo, zeta_upper=zhi, previous_zeta_lower=previous_lo,
        previous_zeta_upper=previous_hi, original_far_weight=weight, kappa=kappa)
    for prefix in ('fixed', 'energy_endpoint'):
        k, fee = (rational(plane, prefix+'_'+key) for key in
                  ('reciprocal_coefficient', 'comparison_fee'))
        values, bounds = plane_floor(rational(plane, 'alpha0'),
            ball(plane, 'timelike_energy'), ball(plane, 'rho_timelike_pairing'),
            ball(plane, 'rho_reference_energy'), ball(plane, 'coordinate_ratio'),
            ball(plane, 'zero_coordinate_energy'), ball(plane, prefix+'_transported_theta'),
            zlo, zhi, rational(plane, 'gamma'), k, fee)
        scalars.update({prefix+'_'+key: value for key, value in values.items()})
        exact.update({prefix+'_'+key: value for key, value in bounds.items()})
        previous = rational(plane, prefix+'_retained_all_alpha_floor_lower')
        exact[prefix+'_previous_floor_lower'] = previous
        exact[prefix+'_retained_floor_lower'] = max(previous, bounds.get('all_alpha_floor_lower', previous))
    scalars.update({key: arb(value) for key, value in exact.items()})
    require_finite(list(scalars.values())+list(left_action.entries())+list(polynomial.coeffs()))
    result = dict(schema='ACTUAL_C27_CACHED_SOURCE_RESPONSE_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={p.name: guards[p] for p in sorted(invoked)},
        inherited_dependency_sha256={name: digest for name, digest in sorted(inherited.items())
                                    if HERE/name not in invoked},
        source_integer_sha256=mu['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=ctx.prec, finite_first_mode=M+1, finite_last_mode=2*M,
        analytic_first_mode=2*M+1, normalization_scale=M, interpolation_order=ORDER,
        interpolation_modes=modes, interpolation_nodes=[[str(x.p), str(x.q)] for x in nodes],
        interpolated_source_coefficients=[envelope(v, BITS) for v in polynomial.coeffs()],
        recovered_approximate_symbols=[envelope(v, BITS) for v in symbols],
        recovered_symbols_are_actual=False, inherited_symbol_error_paid=True,
        physical_left_kernel_action=matrix_out(left_action, BITS),
        scalar_bounds={key: envelope(value, BITS) for key, value in scalars.items()},
        exact_rationals={key: [str(v.p), str(v.q)] for key, v in exact.items()},
        signed_response_strictly_narrower=bool(zlo > previous_lo or zhi < previous_hi),
        fixed_paid_recipe_excluded=bool(exact['fixed_retained_floor_lower'] > 0),
        energy_cap_only_family_excluded=bool(exact['energy_endpoint_retained_floor_lower'] > 0),
        new_cached_symbol_recovery_blocks=1, new_source_interpolations=1,
        new_finite_signed_functionals=1, new_cached_tail_applications=1,
        new_physical_left_kernel_actions=1, new_short_signed_mixed_contractions=1,
        new_paid_plane_applications=2, new_source_sine_phase_moment_inverse_or_solve_suppliers=0,
        old_producers_or_fixtures_replayed=0, parameter_source_precision_cutoff_scan=False,
        actual_center_sign_claim=False, actual_full_border_sign_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.monotonic()-started,
        conditions=['Same real source, original finite Cauchy and approximate old jet, common exact features, physical kernel/lift/reference and weights are inherited, never inferred by interval overlap.',
            'Cauchy support0<n<=M and source l1 upper hold. Lagrange coefficients retain all input uncertainty; remainder holds only on[0,1/2].',
            'Recovered finite symbols enclose the original approximate symbol only; complete global symbol error, tail source error, Abel error, separate physical mixed/forward and smooth errors are paid.',
            'Cached K_q=M^(q-1) sum from2M+1, and phase cap, are reused without any new elementary or moment supplier. Masked physical kernel is not assumed PSD.',
            'Both paid recipes keep their original same actual-row positive-energy upper via saved transported theta. Negative lower is inconclusive; positive floor excludes only the specified paid recipe/family, not an actual smaller-fee center.',
            'No generic theorem, novelty or Lean claim. Actual-center/full-border/good-block/cofinal positivity, Robin and RH remain unresolved.'])
    verify_guards(guards)
    write_packet(out, result)
    print('cached same-source response supplied; seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Require a fresh .zlib output; never overwrite accepted data.')
    ctx.prec = 1024
    consume(args.out)


if __name__ == '__main__':
    main()
