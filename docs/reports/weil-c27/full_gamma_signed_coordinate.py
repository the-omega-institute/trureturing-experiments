"""Supply a raw signed coordinate from the SAME accepted cached source.

Reuse saved symbols, source polynomial, physical left-kernel row and full
physical G00 cap. Only the first lift-column functional and its errors
are new; no old recovery, interpolation, supplier or plane floor runs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, arb_poly, ctx, fmpq
from diagonal_floor_witness import interval_ball
from full_gamma_cached_source_response import finite_pairing, cached_tail
from full_gamma_coherent_split_capacity import freeze_upper
from full_gamma_negative_index_joint import freeze_lower
from full_gamma_signed_mu import require_finite, verify_guards, write_packet
from full_gamma_zero_channel_jet_refine import scaled_forward
from joint_gamma_boundary import matrix_in
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS, M, R = 640, 32768, 56
PINS = {
    'full-gamma-cached-source-response.json.zlib':
        '2f6fa2f46f91dfa0449d5bec2523cdbd91a0e7587284f8e417ccb5bd2068ee03',
    'full-gamma-signed-mu-direct-sine.json.cauchy.zlib':
        'dbfedf8f55f74b4fa1659ff02e57559fea52e3250379c8f99ed1807ee91245f9',
    'full-gamma-signed-mu-direct-sine.json.zlib':
        '6475e5f63dba423b2f1ee7e488581f15ef8d53de83a22ea064ac2e393b124d55',
    'tail-moments.json.zlib':
        'fc91b487a48f82cb276c0823f28c4a38432f12f822d6e0ed2644a67da7c22ad1',
    'full-gamma-lift-projection.json.zlib':
        'a41813b27cd4567af326d42ec77e3b8c8942474b80bdf5833f7bc05aff52f873',
    'full-gamma-signed-correction.json.zlib':
        'd527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd',
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
    'full-gamma-hankel-tail-energy.json.zlib':
        'd09a5ab254d164cfd40982b31e1a155e31d6288b705d916bf248236a0d527606',
    'full-gamma-relaxation-row-witness.json.zlib':
        '1879f4b0653e762ca0d7bb1e1933a59afc79fb10851d155f1c16c87055872992',
    'full-gamma-coordinate-complement.json.zlib':
        '7f6c5d4dc3175e4cd508aadef54d65588bf4c003fddfd01f1c22f8b1748198bb',
}


def coordinate_errors(aa, sigma, weight, delta, ncap, forward,
                      principal, smooth, mixed, kappa, scale=M):
    """Apply original independent symbol/source/physical errors to NEW h0."""
    vals = (sigma, weight, delta, ncap, forward, principal, smooth, mixed, kappa)
    require_finite(list(aa)+[arb(x) for x in vals])
    if (len(aa) != R or type(scale) is not int or scale < 2
            or not weight > 0 or not kappa > 0
            or any(not v >= 0 for v in (sigma, delta, ncap, forward, principal, smooth, mixed))):
        raise ValueError('Missing common symbol/source/physical error domains.')
    ratio = arb(fmpq(scale, scale+1))
    acap = sum((abs(v).upper()*ratio**(2*k) for k, v in enumerate(aa)), arb(0))
    beta = freeze_upper(2*arb(weight).sqrt()*arb(sigma)*acap/(arb.pi()*arb(scale).sqrt()))
    functional = freeze_upper(arb(kappa)*((arb(ncap)+arb(forward))*arb(delta)
                              +(arb(principal)+arb(delta))*arb(beta)))
    smooth_error = freeze_upper(arb(kappa)*arb(smooth)*arb(ncap))
    forward_error = freeze_upper(arb(kappa)*arb(principal)*arb(forward))
    return dict(symbol_jet_error_upper=beta, source_functional_error_upper=functional,
        mixed_cross_error_upper=mixed, smooth_error_upper=smooth_error,
        physical_forward_error_upper=forward_error,
        nonoscillatory_error_upper=functional+mixed+smooth_error+forward_error)


def consume(out):
    started, source = time.perf_counter(), Path(__file__).resolve()
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
            raise ValueError('Different accepted signed-coordinate input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw))
    cached, cauchy, mu, seed, lift, signed, old, hankel, relaxed, complement = packets.values()
    inherited = dict(PINS)
    for packet in packets.values():
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256',
                      'cached_producer_source_sha256'):
            for name, digest in packet.get(field, {}).items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting signed-coordinate input: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked signed-coordinate dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_CACHED_SOURCE_RESPONSE_V1', 'ACTUAL_C27_SIGNED_MU_CAUCHY_BLOCK_V1',
            'ACTUAL_C27_SAME_SOURCE_SIGNED_MU_ABEL_V1', 'ACTUAL_WEIL_NORMALIZED_TAIL_MOMENTS_V1',
            'ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1', 'ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1',
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1', 'ACTUAL_C27_SIGNED_HANKEL_TAIL_ENERGY_V1',
            'C27_SAVED_CONSTRAINT_RELAXATION_ROW_V1', 'ACTUAL_C27_COORDINATE_COMPLEMENT_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any((p['source_bits'], p['source_row'], p['source_integer_sha256']) !=
                   (384, 0, cached['source_integer_sha256'])
                   for p in (cached, cauchy, mu, signed, old, hankel, complement))
            or (relaxed['reference_source_bits'], relaxed['reference_source_row'],
                relaxed['reference_source_integer_sha256']) != (384, 0, cached['source_integer_sha256'])
            or (cached['finite_first_mode'], cached['finite_last_mode'],
                cached['analytic_first_mode'], cached['normalization_scale']) != (M+1, 2*M, 2*M+1, M)
            or (cauchy['finite_first_mode'], cauchy['finite_last_mode'],
                cauchy['normalization_scale']) != (M+1, 2*M, M)
            or cached['recovered_symbols_are_actual'] or not cached['inherited_symbol_error_paid']
            or (seed['c'], seed['M']) != (27, M)
            or seed['analytic_grid_bits'] != mu['symbol_metadata_grid_bits']
            or old['remaining_modes'] != lift['remaining_modes'] or len(lift['remaining_modes']) != 18
            or old['input_sha256'] != lift['input_sha256']
            or (lift['N'], lift['M'], lift['last_mode']) != (1280, M, 2048)
            or lift['basis'] != '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'
            or signed['packet_sha256']['structured-joint-2048-target.json.zlib'] != lift['original_target_sha256']
            or any(complement['packet_sha256'][name] != PINS[name]
                   for name in (tuple(PINS)[0], tuple(PINS)[6], tuple(PINS)[8]))
            or not complement['positive_complement_norm_certified']):
        raise ValueError('Different source, raw lift coordinate, form or cached error contract.')
    if (len(cached['recovered_approximate_symbols']) != M
            or len(cauchy['finite_cauchy_source_outputs']) != M
            or not 1 <= len(cached['interpolated_source_coefficients']) <= 13
            or len(seed['normalized_polynomial']) != 64
            or len(lift['far_features_on_lift']) != 112
            or any(len(row) != 18 for row in lift['far_features_on_lift'])
            or len(cached['physical_left_kernel_action']) != 1
            or len(cached['physical_left_kernel_action'][0]) != 112
            or any(str(q) not in mu['real_hurwitz_kernels'] for q in range(2, 188))):
        raise ValueError('Incomplete cached polynomial/symbol/features/kernel row.')
    rational = lambda p, key: fmpq(*map(int, p['exact_rationals'][key]))
    ball = lambda p, key: interval_ball(p['scalar_bounds'][key], BITS)
    upper = lambda p, key: fmpq(int(p['scalar_bounds'][key][1]), 1 << BITS)
    kappa, weight = rational(cached, 'kappa'), rational(cached, 'original_far_weight')
    if (kappa != fmpq(33, 32) or kappa != rational(signed, 'kappa')
            or kappa != rational(mu, 'kappa') or weight != rational(mu, 'original_far_weight')
            or weight != rational(signed, 'original_far_weight') or not weight > 0):
        raise ValueError('Different original weight or raw-coordinate normalization.')
    # First actual new functional: select e0, never apply P or Q to this coordinate.
    h = matrix_in([[row[0]] for row in lift['far_features_on_lift']], BITS)
    features = list(h.entries())
    new_lift = freeze_upper(interval_ball(lift['lift_gram'][0][0], BITS).upper().sqrt())
    old_lift = rational(signed, 'lift_norm_upper')
    ncap2 = upper(old, 'far_G00_upper')
    if not ncap2 > 0:
        raise ArithmeticError('Missing original physical first-column squared-norm cap.')
    ncap = freeze_upper(arb(ncap2).sqrt())
    forward = scaled_forward(rational(signed, 'sharp_forward_error_upper'), old_lift, new_lift)
    mixed_error = freeze_upper(arb(rational(signed, 'mixed_cross_error_upper'))
                               *arb(new_lift)/arb(old_lift))
    exact = coordinate_errors(features[:R],
        interval_ball(seed['uniform_symbol_error'], seed['analytic_grid_bits']).upper(),
        weight, rational(cached, 'source_tail_norm_error_upper'), ncap, forward,
        upper(hankel, 'principal_weighted_norm_upper'),
        upper(hankel, 'smooth_remainder_weighted_norm_upper'), mixed_error, kappa)
    samples = [interval_ball(v, BITS) for v in cauchy['finite_cauchy_source_outputs']]
    symbols = [interval_ball(v, BITS) for v in cached['recovered_approximate_symbols']]
    polynomial = arb_poly([interval_ball(v, BITS) for v in cached['interpolated_source_coefficients']])
    smooth = arb_poly([interval_ball(v, BITS) for v in seed['normalized_polynomial']])
    kernels = {q: interval_ball(v, BITS) for q, v in mu['real_hurwitz_kernels'].items()}
    finite = finite_pairing(samples, symbols, features, M, kappa, weight)
    tail, tail_cap, oscillatory = cached_tail(polynomial, features, smooth, kernels,
        upper(mu, 'inverse_phase_gap_weight_cap'), M, kappa, weight)
    mixed = arb(kappa)*(matrix_in(cached['physical_left_kernel_action'], BITS)*h)[0, 0]
    error = exact['nonoscillatory_error_upper']+oscillatory
    raw_coordinate = finite+tail+mixed+arb(0, arb(error))
    cap = rational(relaxed, 'relaxed_negative_coordinate_cap')
    lo, hi = max(-cap, freeze_lower(raw_coordinate.lower())), min(cap, freeze_upper(raw_coordinate.upper()))
    if cap < 0 or lo > hi:
        raise ArithmeticError('Signed coordinate conflicts with the SAME-source raw cap; stop.')
    offset = ball(complement, 'plane_coordinate_offset')
    cross_lo = max(rational(complement, 'cross_lower'), lo+freeze_lower(offset.lower()))
    cross_hi = min(rational(complement, 'cross_upper'), hi+freeze_upper(offset.upper()))
    norm_lo, norm_up = [rational(complement, 'projected_norm_'+end) for end in ('lower', 'upper')]
    if cross_lo > cross_hi or not 0 < norm_lo <= norm_up:
        raise ArithmeticError('Inconsistent common-source residual cross or norm.')
    distance = max(fmpq(0), cross_lo, -cross_hi)
    residual = distance*distance/norm_up
    exact.update(new_lift_norm_upper=new_lift, physical_coordinate_norm2_upper=ncap2,
        physical_coordinate_norm_upper=ncap, sharp_forward_error_upper=forward,
        oscillatory_error_upper=oscillatory, complete_coordinate_error_upper=error,
        coordinate_lower=lo, coordinate_upper=hi, cross_lower=cross_lo, cross_upper=cross_hi,
        cross_distance_lower=distance, projected_norm_lower=norm_lo, projected_norm_upper=norm_up,
        residual_lower=residual)
    recipes = {}
    for recipe in ('fixed', 'energy_endpoint'):
        same_f = rational(cached, recipe+'_all_alpha_floor_lower')
        previous = rational(complement, recipe+'_retained_lower')
        augmented, retained = same_f+residual, max(previous, same_f+residual)
        exact.update({recipe+'_same_F_floor_lower':same_f, recipe+'_previous_retained_lower':previous,
                      recipe+'_residual_augmented_lower':augmented, recipe+'_retained_lower':retained})
        recipes[recipe] = dict(strictly_improved=bool(retained > previous),
                               positive_paid_lower_certified=bool(retained > 0))
    scalars = dict(finite_coordinate_functional=finite, tail_coordinate_functional=tail,
        tail_oscillatory_coefficient_cap=tail_cap, mixed_coordinate_kernel_cross=mixed,
        actual_raw_coordinate=raw_coordinate)
    exported = {key: envelope(v, BITS) for key, v in scalars.items()}
    if any(v is None or int(v[0]) > int(v[1]) for v in exported.values()):
        raise ArithmeticError('Invalid outward coordinate export; no retry.')
    result = dict(schema='ACTUAL_C27_SIGNED_RAW_COORDINATE_V1', packet_sha256=PINS,
        source_sha256=guards[source], inherited_dependency_sha256={p.name:h for p,h in guards.items() if p!=source},
        source_integer_sha256=cached['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=ctx.prec, scalar_bounds=exported,
        exact_rationals={key:[str(v.p),str(v.q)] for key,v in exact.items()}, recipes=recipes,
        raw_coordinate_separated=bool(lo > 0 or hi < 0), residual_cross_separated=bool(distance > 0),
        positive_residual_lower_certified=bool(residual > 0), new_coordinate_finite_functionals=1,
        new_coordinate_cached_tail_applications=1, new_short_mixed_dot_products=1,
        new_source_recovery_or_interpolation=0, new_physical_left_kernel_actions=0,
        new_quadratic_matrix_inverse_solve_moment_sine_phase_suppliers=0,
        new_plane_floor_calls=0, old_producers_or_fixtures_replayed=0,
        parameter_source_precision_cutoff_scan=False, actual_center_sign_claim=False,
        actual_full_border_sign_claim=False, actual_rh_proof=False, elapsed_seconds=time.perf_counter()-started,
        conditions=[
            'f=kappa*u_F^*W_F*B, B=R_F V; e0 is the SAME18-coordinate lift first column. f[0] is raw, not divided by sqrt(delta). The saved full G00 bounds ||B e0||_W^2 without kappa.',
            'Cached recovered symbols are approximate, not actual. The source polynomial and physical left-kernel action are reused; no old recovery/interpolation/source/kernel-row action runs.',
            'Original forward and mixed Gram allowances are uniform coefficients times the lift cap. Rescaling by new_L/old_L uses that proof, not an empirical error ratio or a matrix-error-to-vector shortcut.',
            'NEW symbol-jet error, inherited direction-independent source-tail interpolation error, new Abel allowance and separate physical mixed/smooth/forward errors are all paid. Masked kernel is not assumed PSD.',
            'The new signed coordinate and retained offset enclose the SAME actual f and finite form. Their sum is f_perp[0]. Positive saved projected norm permits PSD Cauchy; a zero-containing cross gives residual lower0.',
            'Augment only the saved lower of the SAME F before maximum retention. A nonpositive lower is inconclusive; no actual smaller-fee center/full-border/good-block/cofinal positivity, Robin or RH claim. Existing methods only, no generic novelty or Lean claim.',
        ])
    verify_guards(guards)
    digest = write_packet(out, result)
    print('saved signed raw coordinate', digest, 'seconds', result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Require a fresh .zlib output; never overwrite accepted data.')
    ctx.prec = 1024
    consume(args.out)
