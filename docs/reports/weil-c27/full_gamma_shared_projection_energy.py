"""Bound the SAME source energy from two already saved signed projections.

One new physical short-feature quadratic bounds their common representer.
Classical Gram/Schur bounds retain principal/mixed correlations; all long
sources, projections, inverse profiles and signed suppliers are reused.
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
from full_gamma_mixed_norm_reprice import mixed_norm, reprice
from full_gamma_negative_index_joint import energy_endpoint, freeze_lower, rank_floor
from full_gamma_signed_mu import require_finite, verify_guards, write_packet
from joint_gamma_boundary import matrix_in
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
PINS = {
    'full-gamma-negative-index-joint.json.zlib':
        'a7df6338b629e9dd791ac2914c8efc7f3ef3c69c20e562da0aa18067aa40728d',
    'full-gamma-mixed-norm-reprice.json.zlib':
        '66360516946c571c187b5a4eb35aeed1487eef9eff28bee26d5f8fb3211f1a59',
    'full-gamma-signed-correction.json.zlib':
        'd527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd',
    'full-gamma-signed-mu-direct-sine.json.zlib':
        '6475e5f63dba423b2f1ee7e488581f15ef8d53de83a22ea064ac2e393b124d55',
    'full-gamma-signed-jet-interface.json.zlib':
        '154949d8ab3e1f2875de531aa63a9c5a3c3d7d133640a5ce8b1d23c98d67565a',
    'full-gamma-hankel-tail-energy.json.zlib':
        'd09a5ab254d164cfd40982b31e1a155e31d6288b705d916bf248236a0d527606',
    'structured-joint-2048-target.json.zlib':
        'ea4a761f76b59a11b2d04faf80531901a5164f4145c2424a4e79889ae3f72e22',
    'full-gamma-band-surplus-budget.json.zlib':
        '04fb6eb533b8448ade52ad4032faaaf079ad6603f9998d272a36828c57be44b1',
    'full-gamma-prime-potential-refined.json':
        '9bad7a02ad1a69ff53f0b6a253cce1b04603d1877ba5da05811cbbd141f9389c',
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
}


def projected_energy(p2_lower, p2_upper, c2_lower, c2_upper, norm2_upper, s, t):
    """Real Gram/Schur enclosure with a common representer norm UPPER."""
    require_finite([arb(v) for v in
        (p2_lower, p2_upper, c2_lower, c2_upper, norm2_upper, s, t)])
    if (not 0 <= p2_lower <= p2_upper or not 0 <= c2_lower <= c2_upper
            or not norm2_upper >= 0):
        raise ValueError('Inconsistent squared-norm domains.')
    cauchy = freeze_upper((arb(p2_upper)*arb(c2_upper)).sqrt())
    if norm2_upper == 0:
        if not s.contains(0) or not t.contains(0):
            raise ArithmeticError('A zero representer has nonzero saved projections.')
        cross_lower, cross_upper = -cauchy, cauchy
        residual = arb(cauchy)
    else:
        slo = max(fmpq(0), freeze_lower(abs(s).lower()))
        tlo = max(fmpq(0), freeze_lower(abs(t).lower()))
        ru = freeze_upper(arb(p2_upper)-arb(slo)**2/arb(norm2_upper))
        rc = freeze_upper(arb(c2_upper)-arb(tlo)**2/arb(norm2_upper))
        if min(ru, rc) < 0:
            raise ArithmeticError('Saved projections contradict their common norm cap.')
        residual = (arb(ru)*arb(rc)).sqrt()
        product = s*t/arb(norm2_upper)
        cross_upper = min(cauchy, freeze_upper(product+residual))
        cross_lower = max(-cauchy, freeze_lower(product-residual))
    if cross_lower > cross_upper:
        raise ArithmeticError('Reversed common-projection cross enclosure.')
    total_lower = max(fmpq(0), freeze_lower(arb(p2_lower)+arb(c2_lower)+2*arb(cross_lower)))
    total_upper = freeze_upper(arb(p2_upper)+arb(c2_upper)+2*arb(cross_upper))
    if total_upper < 0 or total_lower > total_upper:
        raise ArithmeticError('Inconsistent joint vector-energy enclosure.')
    return dict(principal_representer_pairing=s, mixed_representer_pairing=t,
                projection_residual_product_upper=residual), dict(
                principal_mixed_cross_lower=cross_lower,
                principal_mixed_cross_upper=cross_upper,
                principal_plus_mixed_energy_lower=total_lower,
                principal_plus_mixed_energy_upper=total_upper)


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
            raise ValueError('Different accepted shared-projection input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    joint, priced, correction, mu, jet, hankel, target, band, refined, old = [
        packets[name] for name in PINS]
    inherited = dict(PINS)
    for packet in (joint, priced, correction):
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256'):
            for name, digest in packet[field].items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting common-projection dependency: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked common-projection dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_NEGATIVE_INDEX_JOINT_BOUND_V1',
            'ACTUAL_C27_SAVED_MIXED_NORM_REPRICE_V1',
            'ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1',
            'ACTUAL_C27_SAME_SOURCE_SIGNED_MU_ABEL_V1',
            'ACTUAL_C27_SIGNED_JET_INTERFACE_V1',
            'ACTUAL_C27_SIGNED_HANKEL_TAIL_ENERGY_V1',
            'ACTUAL_C27_STRUCTURED_TARGET_V1',
            'ACTUAL_C27_BAND_SURPLUS_BUDGET_V1',
            'ACTUAL_C27_REFINED_PRIME_POTENTIAL_BUDGET_V1',
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any(p['source_integer_sha256'] != correction['source_integer_sha256']
                   for p in (joint, priced, mu, jet, hankel, band, refined, old))
            or any((p['source_bits'], p['source_row']) != (384, 0)
                   for p in (joint, priced, correction, mu, jet, hankel, band, old))
            or (jet['physical_last_mode'], jet['first_far_mode'], jet['jet_order'],
                jet['feature_order'], jet['original_metadata_grid_bits']) !=
                (2048, 32769, 56, 'A0..A55,B0..B55', 512)
            or (mu['normalization_scale'], mu['first_source_mode'], mu['last_source_mode'],
                mu['physical_residual_last_mode'], mu['jet_order'], mu['feature_order']) !=
                (32768, 1281, 32768, 2048, 56, 'A0..A55,B0..B55')
            or correction['packet_sha256']['structured-joint-2048-target.json.zlib'] !=
                PINS['structured-joint-2048-target.json.zlib']
            or correction['packet_sha256']['full-gamma-signed-mu-direct-sine.json.zlib'] !=
                PINS['full-gamma-signed-mu-direct-sine.json.zlib']):
        raise ValueError('Different SAME source, representer, physical kernel or scalar units.')
    rational = lambda p, key: fmpq(*map(int, p['exact_rationals'][key]))
    ball = lambda p, key: interval_ball(p['scalar_bounds'][key], BITS)
    upper = lambda p, key: fmpq(int(p['scalar_bounds'][key][1]), 1 << BITS)
    lower = lambda p, key: fmpq(int(p['scalar_bounds'][key][0]), 1 << BITS)
    kappa, weight = rational(correction, 'kappa'), rational(correction, 'original_far_weight')
    eps_far = rational(correction, 'original_far_matrix_error')
    kphys = rational(correction, 'physical_operator_norm_sq_upper')
    lift = rational(correction, 'lift_norm_upper')
    if (kappa != fmpq(33, 32) or kappa != rational(priced, 'kappa')
            or weight != rational(hankel, 'original_far_weight')
            or weight != rational(mu, 'original_far_weight')
            or eps_far != fmpq(int(target['metadata']['far_matrix_error'][1]), 1 << 512)
            or kphys != fmpq(int(target['metadata']['trial_operator_norm_sq_upper'][1]), 1 << 512)
            or lift != rational(jet, 'lift_norm_l1_upper')
            or rational(hankel, 'maximum_reflected_weight') != rational(band, 'maximum_weight')
            or not kphys > 0 or not lift >= 0):
        raise ValueError('Different original physical error, lift norm, weight or charge.')
    h = matrix_in(jet['short_jet_coefficients'], BITS)
    kernel = matrix_in(target['kernel'], BITS)
    if (h.nrows(), h.ncols(), kernel.nrows(), kernel.ncols()) != (112, 1, 112, 112):
        raise ValueError('Incomplete common physical representer features.')
    # Same physical certificate as mixed_norm, now applied to t=D V vR.
    physical_source_norm = freeze_upper(arb(kphys).sqrt()*arb(lift))
    quadratic, error, physical_upper, _ = mixed_norm(h, kernel, eps_far, physical_source_norm)
    old_norm = rational(jet, 'selected_direction_norm_upper')
    norm2 = min(physical_upper, old_norm**2)
    norm = freeze_upper(arb(norm2).sqrt())
    s = ball(mu, 'actual_proxy_functional')/arb(kappa) \
        +arb(0, arb(rational(correction, 'proxy_forward_error_upper')/kappa))
    t = ball(correction, 'mixed_signed_kernel_cross')/arb(kappa) \
        +arb(0, arb(rational(correction, 'mixed_cross_error_upper')/kappa))
    p2lo = max(fmpq(0), freeze_lower(arb(weight*lower(hankel, 'hankel_signed_cross')/4)))
    p2up = freeze_upper(arb(upper(hankel, 'principal_weighted_norm_upper'))**2)
    c2lo = max(fmpq(0), freeze_lower(ball(priced, 'mixed_gram_quadratic')
                                     -arb(rational(priced, 'mixed_gram_error_upper'))))
    c2up = rational(priced, 'mixed_energy_upper')
    scalars, exact = projected_energy(p2lo, p2up, c2lo, c2up, norm2, s, t)
    smooth, constant = [upper(hankel, key) for key in
                        ('smooth_remainder_weighted_norm_upper', 'charge_constant_upper')]
    a1 = upper(refined, 'finite_weighted_first_moment')
    b2 = upper(band, 'high_band_weighted_source_norm_squared_upper')
    maximum, hlo, hhi = rational(band, 'maximum_weight'), rational(old, 'h_lower'), rational(old, 'h_upper')
    # Apply the accepted triangle/pricing map to the NEW split (v+c)+r+0.
    vc_norm = freeze_upper(arb(exact['principal_plus_mixed_energy_upper']).sqrt())
    pricing = reprice(vc_norm, smooth, fmpq(0), rational(priced, 'selected_far_energy_upper'),
                      constant, kappa, weight, a1, b2, maximum, hlo, hhi)
    # The accepted endpoint map now uses ||v+c||^2 >= vc_lower;
    # 4*vc_lower/weight is a rescaled norm certificate, not a Hankel cross.
    endpoint = energy_endpoint(4*exact['principal_plus_mixed_energy_lower']/weight,
        weight, smooth, fmpq(0), constant, kappa, a1, b2, maximum, hlo, hhi)
    if endpoint['source_energy_lower'] <= rational(joint, 'source_energy_lower'):
        endpoint = {key: rational(joint, key) for key in (
            'source_energy_lower', 'endpoint_reciprocal_coefficient_lower', 'endpoint_comparison_fee')}
    if endpoint['source_energy_lower'] > pricing['selected_far_energy_upper']:
        raise ArithmeticError('Conflicting actual source-energy bounds; stop without retry.')
    exact.update(physical_representer_gram_error_upper=error,
        physical_representer_energy_upper=physical_upper, selected_representer_energy_upper=norm2,
        selected_representer_norm_upper=norm, physical_representer_source_norm_upper=physical_source_norm,
        principal_energy_lower=p2lo, principal_energy_upper=p2up,
        mixed_energy_lower=c2lo, mixed_energy_upper=c2up, **pricing, **endpoint)
    for prefix, k, fee in (
            ('fixed', 1/pricing['repriced_legal_coefficient'], pricing['source_valid_fee_upper']),
            ('energy_endpoint', endpoint['endpoint_reciprocal_coefficient_lower'], endpoint['endpoint_comparison_fee'])):
        values = rank_floor(rational(joint, 'timelike_energy_lower'),
            rational(joint, 'rho_pairing_absolute_upper'), rational(joint, 'beta_far_pairing_absolute_upper'), k, fee)
        exact.update({prefix+'_'+key: value for key, value in values.items()})
    scalars.update(physical_representer_gram_quadratic=quadratic,
                   **{key: arb(value) for key, value in exact.items()})
    require_finite(list(scalars.values()))
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    if any(pair is None or int(pair[0]) > int(pair[1]) for pair in exported.values()):
        raise ArithmeticError('Invalid shared-projection export; stop without retry.')
    excluded = lambda prefix: bool(exact.get(prefix+'_all_alpha_floor_lower', fmpq(0)) > 0)
    if excluded('energy_endpoint') and not excluded('fixed'):
        raise ArithmeticError('Contradictory fixed/endpoint comparison.')
    result = dict(schema='ACTUAL_C27_SHARED_PROJECTION_ENERGY_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: digest for p, digest in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: digest for p, digest in guards.items() if p not in invoked},
        source_integer_sha256=correction['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=1024, scalar_bounds=exported,
        exact_rationals={key: [str(value.p), str(value.q)] for key, value in exact.items()},
        strict_representer_norm_improvement=bool(norm < old_norm),
        strict_source_energy_improvement=bool(pricing['selected_far_energy_upper'] < rational(priced, 'selected_far_energy_upper')),
        strict_source_energy_lower_improvement=bool(endpoint['source_energy_lower'] > rational(joint, 'source_energy_lower')),
        fixed_repriced_all_alpha_comparison_excluded=excluded('fixed'),
        fixed_budget_energy_family_excluded=excluded('energy_endpoint'),
        new_representer_feature_quadratics=1, new_principal_mixed_signed_suppliers=0,
        new_source_feature_inverse_trace_or_tail_suppliers=0, old_producers_replayed=0,
        parameter_source_precision_cutoff_scan=False, actual_center_sign_claim=False,
        actual_full_border_sign_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.monotonic()-started,
        conditions=joint['conditions']+[
            'The saved h=P0 t encloses the SAME t=D V vR and ||t||<=sqrt(Kphys)*LV. The original physical Gram error supplies its representer squared norm, without an extra far weight or kappa.',
            'The saved mu/kappa plus the separately paid sharp forward allowance encloses <v,b_actual>; saved nu0/kappa plus its physical mixed-cross allowance encloses <c,b_actual>. All rows and pairings are real and use the SAME original weighted Hilbert space.',
            'For N2>=||b_actual||^2, subtracting the rank-one projection Gram divided by N2 leaves a PSD two-vector Gram. Norms and signed intervals are uniformly enclosed, with no assertion of independent endpoint attainability.',
            'The repricing helper uses the new valid source decomposition (v+c)+r+0. The endpoint helper uses a rescaled joint-norm LOWER certificate, not a new Hankel cross or a far-row upper.',
            'The all-alpha lower comparisons reuse the old accepted B,G,M and their fixed far-row energy UPPER; no constants are rebuilt from an energy lower.'],
        scope='One physical representer quadratic and a saved-projection Gram/Schur application with linked source-energy pricing. No new generic theorem, Lean proof, actual-center/full-border/cofinal, Robin or RH sign.')
    verify_guards(guards)
    digest = write_packet(out, result)
    print('saved', out.name, digest, flush=True)
    for key in ('selected_representer_norm_upper', 'principal_mixed_cross_lower',
                'principal_mixed_cross_upper', 'source_energy_lower', 'selected_far_energy_upper',
                'source_valid_fee_upper', 'repriced_legal_coefficient',
                'fixed_all_alpha_floor_lower', 'energy_endpoint_all_alpha_floor_lower'):
        if key in scalars:
            print(key, scalars[key].str(24), flush=True)
    print('fixed excluded', excluded('fixed'), 'energy family excluded', excluded('energy_endpoint'),
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
