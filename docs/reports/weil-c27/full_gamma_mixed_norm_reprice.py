"""Reprice the SAME source using its saved mixed-H features and Gram.

Supply one missing quadratic, then carry its energy through the source
fee, weighted eligibility cap and affine steering comparison together.
Reuse all old sources, features, kernels, signed responses and bounds.
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
from full_gamma_signed_mu import require_finite, verify_guards, write_packet
from full_gamma_steering_capacity import steering_floor
from joint_gamma_boundary import matrix_in
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
PINS = {
    'full-gamma-signed-correction.json.zlib':
        'd527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd',
    'full-gamma-steering-capacity.json.zlib':
        '85a52adf97e96def95aaf5abaa307338522a2bae5d7b40c96692ad0366a444fe',
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
    'full-gamma-coherent-mixed-capacity.json.zlib':
        'b082eaad102c5f7ca38d1c49c19970952288ce630198a4ea76b48b96d8d23b65',
    'full-gamma-hankel-tail-energy.json.zlib':
        'd09a5ab254d164cfd40982b31e1a155e31d6288b705d916bf248236a0d527606',
    'full-gamma-band-surplus-budget.json.zlib':
        '04fb6eb533b8448ade52ad4032faaaf079ad6603f9998d272a36828c57be44b1',
    'full-gamma-prime-potential-refined.json':
        '9bad7a02ad1a69ff53f0b6a253cce1b04603d1877ba5da05811cbbd141f9389c',
    'structured-joint-2048-target.json.zlib':
        'ea4a761f76b59a11b2d04faf80531901a5164f4145c2424a4e79889ae3f72e22',
}


def mixed_norm(left, kernel, eps_far, lambda_l1):
    """Apply the original physical Gram error to the SAME saved source."""
    size = left.nrows()
    if size < 1 or left.ncols() != 1 or (kernel.nrows(), kernel.ncols()) != (size, size):
        raise ValueError('Incompatible saved feature and physical kernel shapes.')
    require_finite(list(left.entries())+list(kernel.entries())+[arb(eps_far), arb(lambda_l1)])
    if not eps_far >= 0 or not lambda_l1 >= 0:
        raise ValueError('Negative physical error or source norm cap.')
    quadratic = (left.transpose()*kernel*left)[0, 0]
    error = freeze_upper(arb(eps_far)*arb(lambda_l1)**2)
    upper = freeze_upper(quadratic.upper()+arb(error))
    if not upper >= 0:
        raise ArithmeticError('Inconsistent physical mixed-energy upper; stop without retry.')
    norm = freeze_upper(arb(upper).sqrt())
    return quadratic, error, upper, norm


def reprice(principal, smooth, mixed, old_energy, constant, kappa,
            weight, a1, b2, maximum, hlo, hhi):
    """Use ONE common energy cap in the fee and original Bessel budget."""
    vals = [arb(v) for v in (principal, smooth, mixed, old_energy, constant, kappa,
                             weight, a1, b2, maximum, hlo, hhi)]
    require_finite(vals)
    if (any(not v >= 0 for v in (principal, smooth, mixed, b2, maximum))
            or not old_energy > 0 or not kappa > 0 or not weight > 0
            or not a1 > 0 or not 0 < hlo <= hhi):
        raise ValueError('Missing original positive weight, energy or budget domains.')
    rebuilt = freeze_upper((arb(principal)+arb(smooth)+arb(mixed))**2)
    energy = min(old_energy, rebuilt)
    fee = constant+kappa*energy
    budget = freeze_upper(arb(a1)+2*(arb(energy)*arb(b2)/arb(weight)).sqrt()
                          +arb(maximum)*arb(energy)/arb(weight))
    cap = hlo/(hlo+budget)
    coefficient = cap/hhi
    if not fee > 0 or not 0 < cap < 1 or not coefficient > 0:
        raise ArithmeticError('No positive source-valid fee or legal coefficient.')
    return dict(rebuilt_far_energy_upper=rebuilt, selected_far_energy_upper=energy,
                source_valid_fee_upper=fee, complete_weighted_budget_upper=budget,
                supplied_parameter_cap=cap, repriced_legal_coefficient=coefficient)


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
            raise ValueError('Different accepted repricing input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    signed, steering, old, coherent, hankel, band, refined, target = [packets[n] for n in PINS]
    inherited = dict(PINS)
    for packet in (signed, steering):
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256'):
            for name, digest in packet[field].items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting SAME-source dependency: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked repricing dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1',
            'ACTUAL_C27_ALL_ALPHA_COMPARISON_CAPACITY_V1',
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1',
            'ACTUAL_C27_COHERENT_MIXED_CAPACITY_V1',
            'ACTUAL_C27_SIGNED_HANKEL_TAIL_ENERGY_V1',
            'ACTUAL_C27_BAND_SURPLUS_BUDGET_V1',
            'ACTUAL_C27_REFINED_PRIME_POTENTIAL_BUDGET_V1',
            'ACTUAL_C27_STRUCTURED_TARGET_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any(p['source_integer_sha256'] != signed['source_integer_sha256']
                   for p in (steering, old, coherent, hankel, band, refined))
            or any((p['source_bits'], p['source_row']) != (384, 0)
                   for p in (signed, steering, old, coherent, hankel, band))
            or (signed['first_correction_mode'], signed['last_correction_mode'],
                signed['normalization_scale'], signed['jet_order'], signed['feature_order'])
                != (1281, 2048, 32768, 56, 'A0..A55,B0..B55')
            or signed['original_metadata_grid_bits'] != 512
            or signed['packet_sha256']['structured-joint-2048-target.json.zlib'] !=
                PINS['structured-joint-2048-target.json.zlib']):
        raise ValueError('Different source, physical support, scalar grid or packet roles.')
    rational = lambda packet, key: fmpq(*map(int, packet['exact_rationals'][key]))
    upper = lambda packet, key: fmpq(int(packet['scalar_bounds'][key][1]), 1 << BITS)
    ball = lambda packet, key: interval_ball(packet['scalar_bounds'][key], BITS)
    old_energy = rational(old, 'actual_source_energy_cap')
    weight, kappa = rational(hankel, 'original_far_weight'), rational(old, 'kappa')
    eps = rational(signed, 'original_far_matrix_error')
    if (old_energy != upper(hankel, 'far_energy_upper')
            or kappa != fmpq(33, 32) or kappa != rational(signed, 'kappa')
            or weight != rational(signed, 'original_far_weight')
            or rational(hankel, 'maximum_reflected_weight') != rational(band, 'maximum_weight')
            or eps != fmpq(int(target['metadata']['far_matrix_error'][1]), 1 << 512)):
        raise ValueError('Different exact source energy, weight, charge or physical error.')
    left = matrix_in(signed['left_midpoint_features'], BITS)
    kernel = matrix_in(target['kernel'], BITS)
    if left.nrows() != 112 or kernel.nrows() != 112:
        raise ValueError('Incomplete saved mixed feature or kernel.')
    quadratic, gram_error, mixed_energy, mixed_upper = mixed_norm(
        left, kernel, eps, rational(signed, 'lambda_l1_upper'))
    priced = reprice(upper(hankel, 'principal_weighted_norm_upper'),
        upper(hankel, 'smooth_remainder_weighted_norm_upper'), mixed_upper, old_energy,
        upper(hankel, 'charge_constant_upper'), kappa, weight,
        upper(refined, 'finite_weighted_first_moment'),
        upper(band, 'high_band_weighted_source_norm_squared_upper'),
        rational(band, 'maximum_weight'), rational(old, 'h_lower'), rational(old, 'h_upper'))
    energy, fee, tnew = [priced[k] for k in
        ('selected_far_energy_upper', 'source_valid_fee_upper', 'repriced_legal_coefficient')]
    alpha, told = rational(old, 'fixed_alpha'), rational(old, 'legal_coefficient')
    # The accepted producers bound these SAME maps homogeneously in energy.
    ratio = energy/old_energy
    anew = freeze_upper(arb(upper(coherent, 'reference_far_energy_upper'))*arb(ratio))
    cnew = freeze_upper(arb(upper(old, 'negative_far_channel_radius_upper'))*arb(ratio).sqrt())
    h11 = ball(old, 'steering_H11_reference')+arb(1/tnew-1/told)
    r0, lam = ball(old, 'selected_reference_form'), ball(old, 'finite_negative_channel')
    delta = arb(ball(old, 'delta_lower').lower()).union(ball(old, 'delta_upper').upper())
    eta = upper(old, 'neumann_eta_upper')
    if not 0 <= eta < 1:
        raise ValueError('Invalid inherited inverse comparison.')
    inverse_error = freeze_upper(arb(eta)/(1-arb(eta))
        *(r0.upper().sqrt()+arb(anew).sqrt())**2)
    finite = arb(alpha)**2/arb(tnew)+arb(fee)+r0-lam**2/delta
    clin = freeze_upper(finite.upper()+arb(anew)+arb(inverse_error))
    theta = ball(signed, 'actual_signed_theta')
    scalars, exact = steering_floor(alpha, tnew, fee, r0, lam, delta, h11,
        ball(old, 'steering_H12_reference'), ball(old, 'steering_reference_rRr').upper(),
        arb(anew), arb(cnew), theta, arb(eta), clin, arb(inverse_error))
    scalars.update(mixed_gram_quadratic=quadratic, mixed_gram_error_upper=arb(gram_error),
        mixed_energy_upper=arb(mixed_energy), mixed_norm_upper=arb(mixed_upper),
        repriced_reference_far_energy_upper=arb(anew),
        repriced_negative_far_radius_upper=arb(cnew),
        repriced_reference_error_upper=arb(inverse_error),
        center_upper_expression_before_good_metric=arb(clin)+2*theta,
        **{k: arb(v) for k, v in priced.items()})
    exact.update(priced, mixed_gram_error_upper=gram_error, mixed_energy_upper=mixed_energy,
        mixed_norm_upper=mixed_upper, original_far_weight=weight, kappa=kappa,
        old_energy_cap=old_energy, old_legal_coefficient=told,
        old_source_fee=rational(old, 'positive_source_valid_fee'),
        reference_far_energy_upper=anew, negative_far_radius_upper=cnew,
        reference_error_upper=inverse_error)
    require_finite(list(scalars.values()))
    exported = {k: envelope(v, BITS) for k, v in scalars.items()}
    if any(v is None or int(v[0]) > int(v[1]) for v in exported.values()):
        raise ArithmeticError('Invalid repricing export; stop without retry.')
    excluded = scalars['all_alpha_comparison_floor_lower'].lower() > 0
    negative = scalars['center_upper_expression_before_good_metric'].upper() < 0
    if excluded and negative:
        raise ArithmeticError('Contradictory comparison enclosures; stop without retry.')
    result = dict(schema='ACTUAL_C27_SAVED_MIXED_NORM_REPRICE_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: digest for p, digest in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: digest for p, digest in guards.items() if p not in invoked},
        source_integer_sha256=signed['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=1024, scalar_bounds=exported,
        exact_rationals={k: [str(v.p), str(v.q)] for k, v in exact.items()},
        strict_energy_improvement=bool(energy < old_energy),
        strict_fee_improvement=bool(fee < rational(old, 'positive_source_valid_fee')),
        strict_coefficient_improvement=bool(tnew > told),
        all_alpha_repriced_comparison_excluded=bool(excluded),
        negative_upper_before_good_metric=bool(negative),
        actual_center_sign_claim=False, actual_full_border_sign_claim=False, actual_rh_proof=False,
        good_block_metric_paid=False, new_mixed_gram_quadratics=1,
        old_producers_replayed=0, new_source_or_feature_reductions=0,
        new_inverse_reference_trace_gram_or_signed_suppliers=0,
        source_parameter_precision_cutoff_scan=False, elapsed_seconds=time.monotonic()-started,
        conditions=signed['conditions']+steering['conditions']+[
            'The saved a=P0*ell and K0 use the SAME Cphys, exact original Wfar and original physical Gram error; ||ell||<=lambda_l1_upper.',
            'The physical Gram kernel need not be PSD. Its one quadratic plus epsfar*lambda_l1_upper^2 bounds the true mixed squared norm.',
            'uF=u0+rsm+Cphys*ell, and the saved principal/smooth bounds apply to that same source. One common energy cap reprices the independent Kconst upper plus kappa*energy and the original Bessel budget.',
            'The SAME rows rho,betaJ,f and retained lift are independent of the legal source coefficient; changing t0 changes only alpha^2/t0. H11 changes by1/tnew-1/told, and H12, theta0, r0 and rr are unchanged.',
            'The accepted reference positive and negative far caps are homogeneous in this SAME source energy cap; the complete inherited inverse comparison remains valid.',
            'An all-real algebraic comparison includes any legal alpha subset with these fixed data. Positive upper expressions and nonpositive lower capacities determine no actual center sign.'],
        scope='One new source-specific mixed quadratic and a linked fee/eligibility/comparison application of existing bounds. Not an all-energy family certificate, new generic theorem, Lean proof, actual-center/full-border/cofinal, Robin or RH result.')
    verify_guards(guards)
    digest = write_packet(out, result)
    print('saved', out.name, digest, flush=True)
    for key in ('mixed_norm_upper', 'selected_far_energy_upper', 'source_valid_fee_upper',
                'repriced_legal_coefficient', 'center_upper_expression_before_good_metric',
                'all_alpha_comparison_floor_lower'):
        print(key, scalars[key].str(24), flush=True)
    print('negative sufficient upper', negative, 'all-alpha family excluded', excluded,
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
