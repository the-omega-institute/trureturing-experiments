"""Pay the missing jet moments on the SAME saved feature direction.

Reuse the saved h^T K0 h and sharp forward allowance. Signed feature
antidiagonals, three retained moments and directional error payments are new;
all signed projections, source and comparison suppliers are reused.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, ctx, fmpq

from diagonal_floor_witness import interval_ball
from full_gamma_coherent_split_capacity import freeze_upper
from full_gamma_mixed_norm_reprice import reprice
from full_gamma_negative_index_joint import energy_endpoint, freeze_lower, rank_floor
from full_gamma_shared_projection_energy import projected_energy
from full_gamma_signed_mu import require_finite, verify_guards, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS, R, M = 640, 56, 32768
PINS = {
    'full-gamma-shared-projection-energy.json.zlib':
        '59374362635382b3b180aff21584a5a1e32d603b19ccaa5db3ea693086b5f89c',
    'full-gamma-signed-jet-interface.json.zlib':
        '154949d8ab3e1f2875de531aa63a9c5a3c3d7d133640a5ce8b1d23c98d67565a',
    'full-gamma-signed-correction.json.zlib':
        'd527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd',
    'far-coupling.json.zlib':
        'ac35edacd7f55de470586e3b8a5915bcaa49b9d2c4c9c2ee073eeaaa8c9c73e7',
    'full-gamma-primal-tail.json':
        'fad634fb25d8d20e83a15d8277569978646246c2a560d0dc0f7a854b5087455a',
    'full-gamma-hankel-tail-energy.json.zlib':
        'd09a5ab254d164cfd40982b31e1a155e31d6288b705d916bf248236a0d527606',
    'full-gamma-band-surplus-budget.json.zlib':
        '04fb6eb533b8448ade52ad4032faaaf079ad6603f9998d272a36828c57be44b1',
    'full-gamma-prime-potential-refined.json':
        '9bad7a02ad1a69ff53f0b6a253cce1b04603d1877ba5da05811cbbd141f9389c',
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
    'full-gamma-negative-index-joint.json.zlib':
        'a7df6338b629e9dd791ac2914c8efc7f3ef3c69c20e562da0aa18067aa40728d',
}


def omitted_moment_bound(features, symbol_upper, diagonals, anchor, edge):
    """Bound the SAME omitted moments using signed feature antidiagonals."""
    order = len(features)//2
    require_finite(list(features)+list(anchor)+[arb(symbol_upper), arb(edge)])
    if (len(features) != 2*order or order < 1 or type(diagonals) is not int
            or not 1 <= diagonals <= 2*order-1 or len(anchor) != 3
            or not symbol_upper >= 1 or not 0 < edge < 1):
        raise ValueError('Invalid omitted-block features, anchor or support edge.')
    u, t, v = anchor
    determinant = u*v-t*t
    if u.upper() < 0 or v.upper() < 0 or determinant.upper() < 0:
        raise ArithmeticError('Saved low moments have no positive Gram realization.')
    aa, bb = features[:order], features[order:]
    box_total, anchor_total, selected = fmpq(0), fmpq(0), fmpq(0)
    for n in range(diagonals, 2*order-1):
        pairs = [(k, n-k) for k in range(order) if 0 <= n-k < order]
        a = sum((aa[k]*aa[l] for k, l in pairs), arb(0))
        g = sum((aa[k]*bb[l] for k, l in pairs), arb(0))
        b = sum((bb[k]*bb[l] for k, l in pairs), arb(0))
        ap, bp = [max(fmpq(0), freeze_upper(x.upper())) for x in (a, b)]
        gp = freeze_upper(abs(g).upper())
        box = freeze_upper(arb(symbol_upper)**2*arb(ap)/(2*n+1)
            +2*arb(symbol_upper)*arb(gp)/(2*n+2)+arb(bp)/(2*n+3))
        trace_upper = freeze_upper((a*u-2*g*t+b*v).upper())
        delta_lower = freeze_lower(((a*b-g*g)*determinant).lower())
        negative_delta = max(fmpq(0), -delta_lower)
        positive_trace = max(fmpq(0), freeze_upper((arb(trace_upper)
            +(arb(trace_upper)**2+4*arb(negative_delta)).sqrt())/2))
        anchored = freeze_upper(arb(edge)**(2*(n-diagonals+1))*arb(positive_trace))
        box_total += box
        anchor_total += anchored
        selected += min(box, anchored)
    return dict(grouped_box_moment_upper=box_total,
                retained_anchor_moment_upper=anchor_total,
                selected_omitted_moment_upper=selected)


def directional_norm(features, quadratic_upper, symbol_upper, alpha_upper,
                     alpha_radius, weight, moment_error, moment_trace, forward,
                     diagonals, anchor, edge):
    """Enclose the full jet norm from its saved masked quadratic."""
    order = len(features)//2
    values = [arb(v) for v in (quadratic_upper, symbol_upper, alpha_upper,
                               alpha_radius, weight, moment_error, moment_trace, forward)]
    require_finite(list(features)+values)
    if (len(features) != 2*order or order < 1 or type(diagonals) is not int
            or not 1 <= diagonals <= 2*order-1 or not symbol_upper >= 1
            or not alpha_upper > 0 or not weight > 0
            or any(not v >= 0 for v in (alpha_radius, moment_error, moment_trace, forward))):
        raise ValueError('Invalid saved feature, moment, prefactor or forward domains.')
    caps = [freeze_upper(abs(value).upper()) for value in features]
    norm2 = sum((v*v for v in caps), fmpq(0))
    moment_bounds = omitted_moment_bound(features, symbol_upper, diagonals, anchor, edge)
    masking = freeze_upper(arb(weight)*arb(alpha_upper)
                          *arb(moment_bounds['selected_omitted_moment_upper']))
    rounding = freeze_upper(arb(weight)*arb(norm2)*(
        arb(alpha_upper)*arb(moment_error)
        +arb(alpha_radius)*(arb(moment_trace)+arb(moment_error))))
    jet2 = freeze_upper(arb(quadratic_upper)+arb(masking)+arb(rounding))
    if not jet2 >= 0:
        raise ArithmeticError('Inconsistent directional full-jet energy upper.')
    physical2 = freeze_upper((arb(jet2).sqrt()+arb(forward))**2)
    return dict(**moment_bounds, feature_absolute_norm_squared_upper=norm2,
                directional_masking_error_upper=masking,
                directional_rounding_error_upper=rounding,
                full_jet_energy_upper=jet2,
                rebuilt_representer_energy_upper=physical2)


def consume(out):
    started, source = time.monotonic(), Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name and Path(name).resolve().parent == HERE:
            invoked.add(Path(name).resolve())
    guards = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    packets = {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted directional-jet input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    shared, jet, correction, far, primal, hankel, band, refined, old, joint = [
        packets[name] for name in PINS]
    inherited = dict(PINS)
    for packet in (shared, jet, correction):
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256'):
            for name, digest in packet[field].items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting directional-jet dependency: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked directional-jet dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_SHARED_PROJECTION_ENERGY_V1', 'ACTUAL_C27_SIGNED_JET_INTERFACE_V1',
            'ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1', 'ACTUAL_C27_FAR_COUPLING_GRAM_V1',
            'ACTUAL_C27_FULL_GAMMA_PRIMAL_TAIL_V1', 'ACTUAL_C27_SIGNED_HANKEL_TAIL_ENERGY_V1',
            'ACTUAL_C27_BAND_SURPLUS_BUDGET_V1', 'ACTUAL_C27_REFINED_PRIME_POTENTIAL_BUDGET_V1',
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1', 'ACTUAL_C27_NEGATIVE_INDEX_JOINT_BOUND_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any(p['source_integer_sha256'] != shared['source_integer_sha256']
                   for p in (jet, correction, hankel, band, refined, old, joint))
            or any((p['source_bits'], p['source_row']) != (384, 0)
                   for p in (shared, jet, correction, hankel, band, old, joint))
            or (far['c'], far['N'], far['M'], far['jet_order'], far['moment_diagonals'],
                far['bound_grid_bits']) != (27, 1280, M, R, R, 628)
            or (jet['physical_last_mode'], jet['first_far_mode'], jet['jet_order'],
                jet['feature_order'], jet['original_metadata_grid_bits']) !=
                (2048, M+1, R, 'A0..A55,B0..B55', 512)
            or jet['input_sha256']['far-coupling.json.zlib'] != PINS['far-coupling.json.zlib']
            or primal['metadata']['input_sha256'] != jet['input_sha256']
            or shared['packet_sha256']['full-gamma-signed-correction.json.zlib'] !=
                PINS['full-gamma-signed-correction.json.zlib']):
        raise ValueError('Different SAME source, physical jet, moment mask or units.')
    rational = lambda p, key: fmpq(*map(int, p['exact_rationals'][key]))
    upper = lambda p, key: fmpq(int(p['scalar_bounds'][key][1]), 1 << BITS)
    ball = lambda p, key: interval_ball(p['scalar_bounds'][key], BITS)
    width_ok = lambda pair: (len(pair) == 2 and
                             0 <= int(pair[1])-int(pair[0]) <= 2)
    if (len(far['moment_pairs']) != R or len(far['zeta']) != R
            or {row['U_power'] for row in far['moment_pairs']} != set(range(2, 2*R+1, 2))
            or any(row['T_power'] != row['U_power']+1 for row in far['moment_pairs'])
            or {row['p'] for row in far['zeta']} != set(range(4, 2*R+3, 2))
            or any(not width_ok(row[key]) for row in far['moment_pairs'] for key in ('U', 'T'))
            or any(not width_ok(row['interval']) for row in far['zeta'])
            or not width_ok(far['alpha_pair'])):
        raise ValueError('Missing original retained moments or midpoint-radius contract.')
    weight, kappa = rational(correction, 'original_far_weight'), rational(correction, 'kappa')
    hstar = rational(correction, 'moment_gram_norm_upper')
    if (weight != fmpq(int(far['D_numerator']), 1 << far['D_grid_bits'])
            or weight != rational(hankel, 'original_far_weight') or kappa != fmpq(33, 32)
            or hstar != fmpq(int(far['upper_bounds']['moment_gram_norm']), 1 << 628)
            or rational(hankel, 'maximum_reflected_weight') != rational(band, 'maximum_weight')):
        raise ValueError('Different original weight, full moment norm or pricing units.')
    alpha_upper = fmpq(int(far['alpha_pair'][1]), 1 << BITS)
    alpha_radius = fmpq(int(far['alpha_pair'][1])-int(far['alpha_pair'][0]), 1 << (BITS+1))
    moment_error = fmpq(2*R, 1 << BITS)
    # STOP metadata-only review before this first actual feature conversion.
    if len(jet['short_jet_coefficients']) != 2*R or any(
            len(row) != 1 for row in jet['short_jet_coefficients']):
        raise ValueError('Incomplete saved real feature direction.')
    anchor_row = next(row for row in far['moment_pairs'] if row['U_power'] == 2*R)
    anchor_zeta = next(row for row in far['zeta'] if row['p'] == 2*R+2)
    anchor = [interval_ball(pair, BITS) for pair in
              (anchor_row['U'], anchor_row['T'], anchor_zeta['interval'])]
    edge = fmpq(M, M+1)
    features = [interval_ball(row[0], BITS) for row in jet['short_jet_coefficients']]
    symbol = freeze_upper(interval_ball(primal['bounds']['actual_symbol_budget'], BITS).upper())
    exact = directional_norm(features, upper(shared, 'physical_representer_gram_quadratic'),
        symbol, alpha_upper, alpha_radius, weight, moment_error, hstar,
        rational(correction, 'sharp_forward_error_upper'), R, anchor, edge)
    old_norm2 = rational(shared, 'selected_representer_energy_upper')
    norm2 = min(old_norm2, exact['rebuilt_representer_energy_upper'])
    improved = norm2 < old_norm2
    exact.update(selected_representer_energy_upper=norm2,
                 selected_representer_norm_upper=(freeze_upper(arb(norm2).sqrt()) if improved
                     else rational(shared, 'selected_representer_norm_upper')),
                 original_far_weight=weight, kappa=kappa, moment_midpoint_error_upper=moment_error,
                 alpha_upper=alpha_upper, alpha_midpoint_radius=alpha_radius)
    counters = dict(new_projection_applications=0, new_pricing_applications=0, new_rank_floor_applications=0)
    # An ineffective new cap reuses old outputs, without replaying their maps.
    if improved:
        s, z = ball(shared, 'principal_representer_pairing'), ball(shared, 'mixed_representer_pairing')
        _, joint_energy = projected_energy(
            rational(shared, 'principal_energy_lower'), rational(shared, 'principal_energy_upper'),
            rational(shared, 'mixed_energy_lower'), rational(shared, 'mixed_energy_upper'), norm2, s, z)
        exact.update(joint_energy)
        smooth, constant = [upper(hankel, key) for key in
                            ('smooth_remainder_weighted_norm_upper', 'charge_constant_upper')]
        args = (constant, kappa, weight, upper(refined, 'finite_weighted_first_moment'),
                upper(band, 'high_band_weighted_source_norm_squared_upper'),
                rational(band, 'maximum_weight'), rational(old, 'h_lower'), rational(old, 'h_upper'))
        pricing = reprice(freeze_upper(arb(joint_energy['principal_plus_mixed_energy_upper']).sqrt()),
            smooth, fmpq(0), rational(shared, 'selected_far_energy_upper'), *args)
        endpoint = energy_endpoint(4*joint_energy['principal_plus_mixed_energy_lower']/weight,
            weight, smooth, fmpq(0), args[0], args[1], args[3], args[4], args[5], args[6], args[7])
        if endpoint['source_energy_lower'] <= rational(shared, 'source_energy_lower'):
            endpoint = {key: rational(shared, key) for key in (
                'source_energy_lower', 'endpoint_reciprocal_coefficient_lower', 'endpoint_comparison_fee')}
        exact.update(**pricing, **endpoint)
        counters.update(new_projection_applications=1, new_pricing_applications=1)
    else:
        for key in ('principal_mixed_cross_lower', 'principal_mixed_cross_upper',
                    'principal_plus_mixed_energy_lower', 'principal_plus_mixed_energy_upper',
                    'rebuilt_far_energy_upper', 'selected_far_energy_upper', 'source_valid_fee_upper',
                    'complete_weighted_budget_upper', 'supplied_parameter_cap', 'repriced_legal_coefficient',
                    'source_energy_lower', 'endpoint_reciprocal_coefficient_lower', 'endpoint_comparison_fee'):
            exact[key] = rational(shared, key)
    if exact['source_energy_lower'] > exact['selected_far_energy_upper']:
        raise ArithmeticError('Conflicting SAME-source energy bracket; stop without retry.')
    for prefix, k, fee, oldk, oldfee in (
            ('fixed', 1/exact['repriced_legal_coefficient'], exact['source_valid_fee_upper'],
             1/rational(shared, 'repriced_legal_coefficient'), rational(shared, 'source_valid_fee_upper')),
            ('energy_endpoint', exact['endpoint_reciprocal_coefficient_lower'], exact['endpoint_comparison_fee'],
             rational(shared, 'endpoint_reciprocal_coefficient_lower'), rational(shared, 'endpoint_comparison_fee'))):
        if k == oldk and fee == oldfee:
            values = {key: rational(shared, prefix+'_'+key) for key in
                      ('rank_denominator', 'all_alpha_floor_lower') if prefix+'_'+key in shared['exact_rationals']}
        else:
            values = rank_floor(rational(joint, 'timelike_energy_lower'),
                rational(joint, 'rho_pairing_absolute_upper'), rational(joint, 'beta_far_pairing_absolute_upper'), k, fee)
            counters['new_rank_floor_applications'] += 1
        exact.update({prefix+'_'+key: value for key, value in values.items()})
    scalars = {key: arb(value) for key, value in exact.items()}
    require_finite(list(scalars.values()))
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    if any(pair is None or int(pair[0]) > int(pair[1]) for pair in exported.values()):
        raise ArithmeticError('Invalid directional-jet export; stop without retry.')
    excluded = lambda prefix: bool(exact.get(prefix+'_all_alpha_floor_lower', fmpq(0)) > 0)
    if excluded('energy_endpoint') and not excluded('fixed'):
        raise ArithmeticError('Contradictory fixed/endpoint comparison.')
    result = dict(schema='ACTUAL_C27_DIRECTIONAL_JET_NORM_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: digest for p, digest in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: digest for p, digest in guards.items() if p not in invoked},
        source_integer_sha256=shared['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=1024, scalar_bounds=exported,
        exact_rationals={key: [str(value.p), str(value.q)] for key, value in exact.items()},
        strict_representer_norm_improvement=improved,
        strict_source_energy_upper_improvement=bool(exact['selected_far_energy_upper'] < rational(shared, 'selected_far_energy_upper')),
        strict_source_energy_lower_improvement=bool(exact['source_energy_lower'] > rational(shared, 'source_energy_lower')),
        fixed_repriced_all_alpha_comparison_excluded=excluded('fixed'),
        fixed_budget_energy_family_excluded=excluded('energy_endpoint'),
        new_omitted_block_contractions=1, new_feature_absolute_norms=1,
        new_signed_feature_antidiagonal_families=3, retained_anchor_scalar_reads=3,
        new_masked_feature_quadratics=0, new_signed_or_moment_suppliers=0,
        new_source_feature_inverse_trace_or_tail_suppliers=0, old_producers_replayed=0,
        parameter_source_precision_cutoff_scan=False, **counters,
        actual_center_sign_claim=False, actual_full_border_sign_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.monotonic()-started,
        conditions=shared['conditions']+[
            'SAME full jet j_h has Gram alpha*wF*H and the saved h=P0t. The saved Qupper encloses h^T(alpha0*wF*H0)h; no quadratic is repeated.',
            'The missing block k+l>=56 is grouped by signed feature antidiagonals before taking bounds. Positive AA/BB coefficients and absolute AB coefficients give the box bound, with both AB blocks counted.',
            'The three retained true moments U112,T113,V114 enclose a PSD Gram Sstar. Actual support x=M/m<=M/(M+1) gives 0<=Sn<=edge^(2(n-55))*Sstar. The two-by-two positive-part trace formula gives a second upper; the per-antidiagonal minimum is uniform on the SAME realization. Rounded centers need not be PSD.',
            'The original exact midpoint H0 and alpha0 use640-grid intervals. Moment spectral error <=112*2^-640 and trace(H)<=Hstar; masking reduces Frobenius norm, hence ||H0||<=Hstar+moment_error justify the separate directional rounding debit.',
            'The existing sharp forward norm allowance pays b_actual-j_h. The complete physical Gram error times ||t||^2 is replaced by these direct jet payments, not silently dropped.',
            'Only a strictly improved representer cap is sent to the already accepted projection/pricing maps. Unchanged fixed or endpoint maps are reused from their saved outputs. B,G,M and their far-row energy UPPER stay fixed.'],
        scope='Directional application of existing exact jet Gram, midpoint and moment bounds. No new generic theorem, Lean, novelty, actual-center/full-border/cofinal, Robin or RH sign.')
    verify_guards(guards)
    digest = write_packet(out, result)
    print('saved', out.name, digest, flush=True)
    for key in ('directional_masking_error_upper', 'directional_rounding_error_upper',
                'selected_representer_norm_upper', 'principal_mixed_cross_upper',
                'source_energy_lower', 'selected_far_energy_upper', 'source_valid_fee_upper',
                'repriced_legal_coefficient', 'fixed_all_alpha_floor_lower',
                'energy_endpoint_all_alpha_floor_lower'):
        if key in scalars:
            print(key, scalars[key].str(24), flush=True)
    print('norm improved', improved, 'fixed excluded', excluded('fixed'),
          'energy family excluded', excluded('energy_endpoint'),
          'seconds', time.monotonic()-started, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output. Fixed1024/640; no scan or retry.')
    ctx.prec = 1024
    consume(args.out)


if __name__ == '__main__':
    main()
