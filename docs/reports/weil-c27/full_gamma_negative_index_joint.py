"""Use a saved timelike row to bound the SAME repriced comparison.

Apply the classical negative-index-one inequality and scalar square
completion. Consume saved scalars only; no derivative or Gram supplier.
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
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
PINS = {
    'full-gamma-mixed-norm-reprice.json.zlib':
        '66360516946c571c187b5a4eb35aeed1487eef9eff28bee26d5f8fb3211f1a59',
    'full-gamma-steering-capacity.json.zlib':
        '85a52adf97e96def95aaf5abaa307338522a2bae5d7b40c96692ad0366a444fe',
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
    'full-gamma-signed-correction.json.zlib':
        'd527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd',
    'full-gamma-hankel-tail-energy.json.zlib':
        'd09a5ab254d164cfd40982b31e1a155e31d6288b705d916bf248236a0d527606',
    'full-gamma-band-surplus-budget.json.zlib':
        '04fb6eb533b8448ade52ad4032faaaf079ad6603f9998d272a36828c57be44b1',
    'full-gamma-prime-potential-refined.json':
        '9bad7a02ad1a69ff53f0b6a253cce1b04603d1877ba5da05811cbbd141f9389c',
}


def freeze_lower(value):
    pair = envelope(arb(value), BITS)
    if pair is None or int(pair[0]) > int(pair[1]):
        raise ArithmeticError('Nonfinite lower export; stop without retry.')
    return fmpq(int(pair[0]), 1 << BITS)


def joint_constants(alpha, told, fee_old, qold, h11, h12, r0, rr, theta, eta, afar):
    """Enclose b, rho Qeps y0* and (beta+f) Qeps y0* jointly."""
    require_finite([arb(v) for v in
        (alpha, told, fee_old, qold, h11, h12, r0, rr, theta, eta, afar)])
    if (not told > 0 or not 0 <= eta < 1
            or any(not v >= 0 for v in (r0, rr, afar))):
        raise ValueError('Missing inverse or positive-seminorm domains.')
    eps = eta/(1+eta)
    d0 = h11*arb(alpha)+h12
    b = arb(fee_old)+arb(alpha)**2/arb(told)-qold+arb(eps)*r0
    positive_cross = (arb(rr)*arb(r0.upper())).sqrt().upper()
    gref = d0-arb(alpha)/arb(told)
    cref = qold-arb(fee_old)+theta-arb(alpha)*d0
    g = freeze_upper(abs(gref).upper()+arb(eps)*positive_cross)
    m = freeze_upper(abs(cref).upper()+arb(eps)*(r0.upper()
        +abs(arb(alpha)).upper()*positive_cross+(arb(afar)*r0.upper()).sqrt().upper()))
    lower = freeze_lower(b)
    if not lower > 0:
        raise ArithmeticError('No saved timelike reference; stop without retry.')
    return dict(timelike_energy=b, reference_derivative=d0,
                reference_rho_pairing=gref, reference_beta_far_pairing=cref), dict(
                epsilon=eps, timelike_energy_lower=lower,
                rho_pairing_absolute_upper=g, beta_far_pairing_absolute_upper=m)


def rank_floor(b, g, m, k, fee):
    """Complete k alpha^2+fee-(g|alpha|+m)^2/b over all real alpha."""
    if not b > 0 or not k > 0 or min(g, m) < 0:
        raise ValueError('Missing rank-one scalar domains.')
    denominator = b*k-g*g
    values = dict(rank_denominator=denominator)
    if denominator > 0:
        values['all_alpha_floor_lower'] = fee-k*m*m/denominator
    return values


def energy_endpoint(cross_lower, weight, smooth_upper, mixed_upper,
                    constant, kappa, a1, b2, maximum, hlo, hhi):
    """A source-energy LOWER endpoint; it never bounds the fixed far row."""
    require_finite([arb(v) for v in (cross_lower, weight, smooth_upper, mixed_upper,
        constant, kappa, a1, b2, maximum, hlo, hhi)])
    if (not weight > 0 or not kappa > 0 or not a1 > 0
            or not 0 < hlo <= hhi
            or min(smooth_upper, mixed_upper, b2, maximum) < 0):
        raise ValueError('Missing original weighted-energy or Bessel domains.')
    principal_lower = (arb(max(weight*cross_lower/4, fmpq(0)))).sqrt().lower()
    norm_lower = max(freeze_lower(principal_lower-arb(smooth_upper)-arb(mixed_upper)), fmpq(0))
    energy_lower = freeze_lower(arb(norm_lower)**2)
    # Unrounded polynomial budget is <= every upward-rounded recipe budget.
    budget = arb(a1)+2*(arb(energy_lower)*arb(b2)/arb(weight)).sqrt() \
        +arb(maximum)*arb(energy_lower)/arb(weight)
    k_lower = freeze_lower(arb(hhi)*(1+budget/arb(hlo)))
    return dict(source_energy_lower=energy_lower, endpoint_reciprocal_coefficient_lower=k_lower,
                endpoint_comparison_fee=constant+kappa*energy_lower)


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
            raise ValueError('Different accepted joint-bound input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    priced, steering, old, signed, hankel, band, refined = [packets[n] for n in PINS]
    inherited = dict(PINS)
    for packet in (priced, steering):
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256'):
            for name, digest in packet[field].items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting joint-bound dependency: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked joint-bound dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_SAVED_MIXED_NORM_REPRICE_V1',
            'ACTUAL_C27_ALL_ALPHA_COMPARISON_CAPACITY_V1',
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1',
            'ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1',
            'ACTUAL_C27_SIGNED_HANKEL_TAIL_ENERGY_V1',
            'ACTUAL_C27_BAND_SURPLUS_BUDGET_V1',
            'ACTUAL_C27_REFINED_PRIME_POTENTIAL_BUDGET_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any(p['source_integer_sha256'] != signed['source_integer_sha256']
                   for p in packets.values())
            or any((p['source_bits'], p['source_row']) != (384, 0)
                   for p in (priced, steering, old, signed, hankel, band))):
        raise ValueError('Different SAME-source scalar roles or units.')
    rational = lambda p, key: fmpq(*map(int, p['exact_rationals'][key]))
    ball = lambda p, key: interval_ball(p['scalar_bounds'][key], BITS)
    upper = lambda p, key: fmpq(int(p['scalar_bounds'][key][1]), 1 << BITS)
    lower = lambda p, key: fmpq(int(p['scalar_bounds'][key][0]), 1 << BITS)
    alpha, told, fee_old = [rational(old, key) for key in
        ('fixed_alpha', 'legal_coefficient', 'positive_source_valid_fee')]
    if (rational(priced, 'alpha0') != alpha
            or rational(priced, 'old_legal_coefficient') != told
            or rational(priced, 'old_source_fee') != fee_old
            or rational(steering, 'alpha0') != alpha
            or rational(steering, 'legal_coefficient') != told
            or rational(steering, 'frozen_positive_fee') != fee_old
            or rational(priced, 'original_far_weight') != rational(hankel, 'original_far_weight')
            or rational(priced, 'kappa') != fmpq(33, 32)
            or rational(hankel, 'maximum_reflected_weight') != rational(band, 'maximum_weight')):
        raise ValueError('Different fixed reference, weight or repricing inputs.')
    scalars, exact = joint_constants(alpha, told, fee_old,
        ball(steering, 'reference_finite_expression'), ball(old, 'steering_H11_reference'),
        ball(old, 'steering_H12_reference'), ball(old, 'selected_reference_form'),
        upper(old, 'steering_reference_rRr'), ball(signed, 'actual_signed_theta'),
        upper(old, 'neumann_eta_upper'), rational(priced, 'reference_far_energy_upper'))
    endpoint = energy_endpoint(lower(hankel, 'hankel_signed_cross'),
        rational(hankel, 'original_far_weight'), upper(hankel, 'smooth_remainder_weighted_norm_upper'),
        rational(priced, 'mixed_norm_upper'), upper(hankel, 'charge_constant_upper'),
        rational(priced, 'kappa'), upper(refined, 'finite_weighted_first_moment'),
        upper(band, 'high_band_weighted_source_norm_squared_upper'),
        rational(band, 'maximum_weight'), rational(old, 'h_lower'), rational(old, 'h_upper'))
    knew = 1/rational(priced, 'repriced_legal_coefficient')
    dnew = rational(priced, 'source_valid_fee_upper')
    if (not 0 <= endpoint['source_energy_lower'] <= rational(priced, 'selected_far_energy_upper')
            or not 0 < endpoint['endpoint_reciprocal_coefficient_lower'] <= knew
            or not endpoint['endpoint_comparison_fee'] <= dnew):
        raise ArithmeticError('Inconsistent joint energy endpoint; stop without retry.')
    exact.update(endpoint, fixed_reciprocal_coefficient=knew, fixed_comparison_fee=dnew,
                 fixed_far_positive_energy_upper=rational(priced, 'reference_far_energy_upper'))
    for prefix, k, fee in (
            ('fixed', knew, dnew),
            ('energy_endpoint', endpoint['endpoint_reciprocal_coefficient_lower'],
             endpoint['endpoint_comparison_fee'])):
        values = rank_floor(exact['timelike_energy_lower'], exact['rho_pairing_absolute_upper'],
                            exact['beta_far_pairing_absolute_upper'], k, fee)
        exact.update({prefix+'_'+key: value for key, value in values.items()})
    scalars.update({key: arb(value) for key, value in exact.items()})
    require_finite(list(scalars.values()))
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    if any(pair is None or int(pair[0]) > int(pair[1]) for pair in exported.values()):
        raise ArithmeticError('Invalid scalar joint-bound export; stop without retry.')
    excluded = lambda prefix: bool(exact.get(prefix+'_all_alpha_floor_lower', fmpq(0)) > 0)
    if excluded('energy_endpoint') and not excluded('fixed'):
        raise ArithmeticError('Contradictory monotone endpoint comparison.')
    result = dict(schema='ACTUAL_C27_NEGATIVE_INDEX_JOINT_BOUND_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: digest for p, digest in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: digest for p, digest in guards.items() if p not in invoked},
        source_integer_sha256=signed['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=1024, scalar_bounds=exported,
        exact_rationals={key: [str(value.p), str(value.q)] for key, value in exact.items()},
        fixed_repriced_all_alpha_comparison_excluded=excluded('fixed'),
        fixed_budget_energy_family_excluded=excluded('energy_endpoint'),
        new_scalar_joint_application=True, new_derivative_or_gram_suppliers=0,
        old_producers_replayed=0, parameter_source_precision_cutoff_scan=False,
        actual_center_sign_claim=False, actual_full_border_sign_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.monotonic()-started,
        conditions=priced['conditions']+[
            'P=T R T* is PSD; Qeps=(1-eps)P-e0e0*/delta <= Qtrue for eps=eta/(1+eta). The saved y0 is timelike for this SAME Qeps.',
            'qold, d0 and theta retain their common reference realization. Rank-one negative index bounds x Qeps x* below by -|x Qeps y0*|^2/b.',
            'The far-row seminorm cap in B,G,M remains frozen at the accepted NEW energy UPPER, even at the artificial energy LOWER endpoint.',
            'The signed Hankel cross lower is part of the SAME weighted principal norm squared, normalized by1/4; the remaining diagonal is nonnegative.',
            'The whole-energy claim, when the endpoint passes, varies only the energy cap in D=Kconst_up+kappa E and the original upward-rounded Bessel budget with fixed first-moment/b2/maximum/h endpoints. It does not cover other fee or budget improvements.',
            'A positive comparison floor excludes this paid recipe, not an actual smaller-fee center, another source/lift, full-border positivity, Robin or RH.'],
        scope='Saved-scalar application of classical negative-index-one/Schur and square-completion bounds. No new generic theorem, Lean proof or actual-center/full-border/cofinal/RH sign.')
    verify_guards(guards)
    digest = write_packet(out, result)
    print('saved', out.name, digest, flush=True)
    for key in ('timelike_energy_lower', 'rho_pairing_absolute_upper',
                'beta_far_pairing_absolute_upper', 'source_energy_lower',
                'fixed_rank_denominator', 'fixed_all_alpha_floor_lower',
                'energy_endpoint_rank_denominator', 'energy_endpoint_all_alpha_floor_lower'):
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
