"""Test one equality row in the saved c27 scalar relaxation.

The constructed row is not the actual far response. Reuse the finite core,
signed-response interval, energy caps and inverse comparison; no physical
source, derivative, Gram, inverse or long-tail supplier is evaluated.
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
ALPHA = fmpq(101, 2000)
PINS = {
    'full-gamma-directional-jet-norm.json.zlib':
        '9ccd6494c0520f9f53efe0b44bcae879566157c7e7ffb698165f9d02c258fd41',
    'full-gamma-mixed-norm-reprice.json.zlib':
        '66360516946c571c187b5a4eb35aeed1487eef9eff28bee26d5f8fb3211f1a59',
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
    'full-gamma-negative-index-joint.json.zlib':
        'a7df6338b629e9dd791ac2914c8efc7f3ef3c69c20e562da0aa18067aa40728d',
    'full-gamma-signed-correction.json.zlib':
        'd527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd',
}


def relaxation_row(alpha, alpha0, theta, delta, lam, rr, rb, r0,
                   gref, eta, k, fee, positive_cap, negative_cap):
    """One fixed-alpha timelike equality row under the original core interface."""
    vals = [arb(x) for x in (alpha, alpha0, theta, delta, lam, rr, rb,
                            r0, gref, eta, k, fee, positive_cap, negative_cap)]
    require_finite(vals)
    if (not delta.lower() > 0 or lam.contains(0)
            or not rr.lower() >= 0 or not r0.lower() >= 0
            or not 0 <= eta < 1 or not k > 0 or not fee > 0
            or min(positive_cap, negative_cap) < 0):
        raise ValueError('Missing finite-form, nonzero coordinate or cap domains.')
    p = arb(alpha0)*rr+rb
    b = lam**2/delta-r0
    if not b.lower() > 0:
        raise ArithmeticError('The saved reference is not certified timelike for Q_R.')
    s = arb(alpha)-arb(alpha0)
    tau = -(arb(theta)+s*gref)/b
    # For the actual finite parameters define f=tau*y0-s*rho.
    # Then f Q_R y0*=theta and alpha*rho+beta+f=(1+tau)*y0.
    positive_energy = tau**2*r0-2*tau*s*p+s**2*rr
    rho0 = delta*(p-gref)/lam
    negative_coordinate = tau*lam-s*rho0
    if positive_energy.upper() < 0:
        raise ArithmeticError('Inconsistent positive reference-form interval.')
    energy_upper = freeze_upper(positive_energy.upper())
    coordinate_upper = freeze_upper(abs(negative_coordinate).upper())
    eps_plus = eta/(1-eta)
    # Q_true <= Q_R+eps_plus*P on this same finite core.
    expression = (arb(k)*arb(alpha)**2+arb(fee)
        +(1+tau)**2*(-b+arb(eps_plus)*r0))
    expression_upper = freeze_upper(expression.upper())
    scalars = dict(reference_timelike_energy=b,
        reference_rho_positive_pairing=p, relaxed_row_coefficient=tau,
        relaxed_row_positive_energy=positive_energy,
        relaxed_row_negative_coordinate=negative_coordinate,
        relaxed_paid_comparison_upper_expression=expression)
    exact = dict(relaxed_positive_energy_upper=energy_upper,
        relaxed_negative_coordinate_absolute_upper=coordinate_upper,
        relaxed_paid_comparison_upper=expression_upper,
        inverse_upper_relative_error=eps_plus)
    flags = dict(relaxed_positive_cap_satisfied=bool(energy_upper <= positive_cap),
        relaxed_negative_cap_satisfied=bool(coordinate_upper <= negative_cap),
        relaxed_comparison_strictly_negative=bool(expression_upper < 0))
    flags['relaxed_negative_row_certified'] = all(flags.values())
    return scalars, exact, flags


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
            raise ValueError('Different accepted relaxation input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw))
    new, priced, old, joint, signed = [packets[name] for name in PINS]
    inherited = dict(PINS)
    for packet in (new, priced, joint):
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256'):
            for name, digest in packet[field].items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting saved relaxation dependency: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked relaxation dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_DIRECTIONAL_JET_NORM_V1',
            'ACTUAL_C27_SAVED_MIXED_NORM_REPRICE_V1',
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1',
            'ACTUAL_C27_NEGATIVE_INDEX_JOINT_BOUND_V1',
            'ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any(p['source_integer_sha256'] != new['source_integer_sha256']
                   for p in packets.values())
            or any((p['source_bits'], p['source_row']) != (384, 0)
                   for p in packets.values())):
        raise ValueError('Different reference source, signed response or scalar units.')
    rational = lambda p, key: fmpq(*map(int, p['exact_rationals'][key]))
    upper = lambda p, key: fmpq(int(p['scalar_bounds'][key][1]), 1 << BITS)
    ball = lambda p, key: interval_ball(p['scalar_bounds'][key], BITS)
    alpha0 = rational(old, 'fixed_alpha')
    if (rational(priced, 'alpha0') != alpha0
            or rational(new, 'original_far_weight') != rational(priced, 'original_far_weight')
            or rational(new, 'kappa') != rational(priced, 'kappa')
            or rational(new, 'kappa') != fmpq(33, 32)
            or joint['packet_sha256']['full-gamma-negative-channel-capacity.json.zlib'] !=
                PINS['full-gamma-negative-channel-capacity.json.zlib']):
        raise ValueError('Different finite frame or homogeneous far-cap interface.')
    old_energy = rational(priced, 'selected_far_energy_upper')
    new_energy = rational(new, 'selected_far_energy_upper')
    if not 0 < new_energy <= old_energy:
        raise ValueError('Missing same-source energy upper improvement.')
    ratio = new_energy/old_energy
    positive_cap = freeze_upper(arb(rational(priced, 'reference_far_energy_upper'))*arb(ratio))
    negative_cap = freeze_upper(arb(rational(priced, 'negative_far_radius_upper'))*arb(ratio).sqrt())
    theta_pair = signed['scalar_bounds']['actual_signed_theta']
    if not int(theta_pair[0]) < int(theta_pair[1]):
        raise ValueError('No interior response value in the saved interval.')
    # This target is chosen inside the interval, not supplied as the actual response.
    target = fmpq(int(theta_pair[0])+int(theta_pair[1]), 1 << (BITS+1))
    k, fee = 1/rational(new, 'repriced_legal_coefficient'), rational(new, 'source_valid_fee_upper')
    eta = upper(old, 'neumann_eta_upper')
    scalars, exact, flags = relaxation_row(ALPHA, alpha0, target,
        -ball(old, 'old_first_pivot'), ball(old, 'finite_negative_channel'),
        ball(old, 'steering_reference_rRr'), ball(old, 'steering_reference_rRb'),
        ball(old, 'selected_reference_form'), ball(joint, 'reference_rho_pairing'),
        eta, k, fee, positive_cap, negative_cap)
    exact.update(fixed_relaxation_alpha=ALPHA, reference_alpha0=alpha0,
        chosen_response_target=target, comparison_reciprocal_coefficient=k,
        comparison_fee=fee, source_energy_upper=new_energy,
        reference_source_energy_upper=old_energy, homogeneous_energy_ratio=ratio,
        relaxed_positive_energy_cap=positive_cap, relaxed_negative_coordinate_cap=negative_cap,
        inverse_error_eta_upper=eta)
    scalars.update({key: arb(value) for key, value in exact.items()})
    require_finite(list(scalars.values()))
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    if any(pair is None or int(pair[0]) > int(pair[1]) for pair in exported.values()):
        raise ArithmeticError('Invalid relaxation-row export; stop without retry.')
    result = dict(schema='C27_SAVED_CONSTRAINT_RELAXATION_ROW_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: digest for p, digest in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: digest for p, digest in guards.items() if p not in invoked},
        reference_source_integer_sha256=new['source_integer_sha256'],
        reference_source_bits=384, reference_source_row=0, grid_bits=BITS, working_bits=1024,
        scalar_bounds=exported,
        exact_rationals={key: [str(value.p), str(value.q)] for key, value in exact.items()},
        **flags, chosen_target_is_actual_signed_response=False,
        physical_source_realization_claim=False, physical_steering_legality_claim=False,
        new_relaxation_row_applications=1, new_homogeneous_cap_rescalings=1,
        new_actual_source_derivative_gram_inverse_or_tail_suppliers=0,
        old_producers_replayed=0, parameter_source_precision_cutoff_scan=False,
        actual_center_sign_claim=False, actual_full_border_sign_claim=False, actual_rh_proof=False,
        reference_conditions=new['conditions'],
        conditions=[
            'P=T R T* and Q_R=P-e0e0*/delta are the SAME finite reference form. Saved rr,rb,r0,lambda and gref enclose their common values; lambda=y0[0] is nonzero and y0=alpha0*rho+beta.',
            'The original two-sided inverse comparison gives Qtrue<=Q_R+eta/(1-eta)*P. This is an upper used on the constructed row, never a lower-to-upper substitution.',
            'The original far positive-energy and negative-coordinate caps are homogeneous in the SAME source-energy UPPER. The accepted rescaling uses the ratio of new and old uppers, not an energy lower.',
            'Choose the one literal alpha101/2000 and an interior theta target. Define s=alpha-alpha0,b=-y0 Q_R y0*,tau=-(theta+s*gref)/b and f_rel=tau*y0-s*rho. This f_rel is an algebraic row in a scalar relaxation, not the actual far response.',
            'Then f_rel Q_R y0*=theta exactly, and alpha*rho+beta+f_rel=(1+tau)*y0. Its positive norm and negative coordinate are explicitly bounded against both saved caps before claiming a negative relaxed comparison.',
            'The chosen theta is inside the measured interval, but need not equal the actual theta. Membership of the chosen alpha in a physical legal subset and realization of f_rel by the original source/coupling are not asserted.',
            'Even a certified negative relaxed comparison is not a negative value for the original source, an actual-center sign, full-border/cofinal positivity, Robin or RH.'
        ],
        scope='One saved-scalar application of classical timelike Schur equality and inverse upper comparison to a specified relaxation. No new generic theorem, Lean proof, actual source or RH claim.',
        elapsed_seconds=time.monotonic()-started)
    verify_guards(guards)
    digest = write_packet(out, result)
    print('saved', out.name, digest, flush=True)
    for key in ('relaxed_row_coefficient', 'relaxed_positive_energy_upper',
                'relaxed_positive_energy_cap', 'relaxed_negative_coordinate_absolute_upper',
                'relaxed_negative_coordinate_cap', 'relaxed_paid_comparison_upper'):
        print(key, scalars[key].str(24), flush=True)
    print('relaxation only', flags, 'seconds', time.monotonic()-started, flush=True)


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
