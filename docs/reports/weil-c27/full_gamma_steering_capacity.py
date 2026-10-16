"""Test all real steering coefficients using only saved scalar certificates.

Keep the original source, legal coefficient, fee, old core and lift.
Price a uniform comparison-expression floor before any new signed
derivative, source, feature, matrix, inverse or tail computation.
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
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
    'full-gamma-coherent-mixed-capacity.json.zlib':
        'b082eaad102c5f7ca38d1c49c19970952288ce630198a4ea76b48b96d8d23b65',
    'full-gamma-signed-correction.json.zlib':
        'd527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd',
}


def steering_floor(alpha, t0, fee, r0, lam, delta, h11, h12, rr,
                   positive_far, cminus, theta, eta, clin, old_error):
    """Apply Cauchy, the saved inverse comparison and scalar completion."""
    values = [arb(v) for v in (alpha, t0, fee, r0, lam, delta, h11, h12, rr,
                               positive_far, cminus, theta, eta, clin, old_error)]
    require_finite(values)
    if (not t0 > 0 or not fee > 0 or not delta.lower() > 0
            or not h11.lower() > 0 or not 0 <= eta < 1
            or any(not v >= 0 for v in (rr, positive_far, cminus, old_error))):
        raise ValueError('Missing positive curvature/coefficient or norm domains.')
    # h11=1/t0+rr_exact-rho0^2/delta, with rr an upper for rr_exact.
    negative_rho = (1/arb(t0)+arb(rr)-h11.lower()).upper()
    if not negative_rho >= 0:
        raise ArithmeticError('Inconsistent saved curvature and positive reference norm.')
    negative_rho_up = freeze_upper(negative_rho)
    fcredit = freeze_upper(arb(cminus)**2/delta.lower())
    phi_positive = (arb(positive_far)*arb(rr)).sqrt().upper()
    phi_negative = (arb(fcredit)*arb(negative_rho_up)).sqrt().upper()
    phi_cap = freeze_upper(phi_positive+phi_negative)
    derivative = h11*arb(alpha)+h12
    slope_cap = freeze_upper(abs(derivative).upper()+arb(phi_cap))
    finite = arb(alpha)**2/arb(t0)+arb(fee)+r0-lam**2/delta
    # Qtrue >= Qref-eps*Qpositive_ref; both comparison directions are saved.
    eps = (arb(eta)/(1+arb(eta))).upper()
    # ||y0+s*rho+f||_positive^2 <= 4*r0_up+4*s^2*rr+2*A.
    constant_loss = freeze_upper(eps*(4*r0.upper()+2*arb(positive_far)))
    curvature_loss = freeze_upper(4*eps*arb(rr))
    curvature = h11.lower()-arb(curvature_loss)
    if not curvature > 0:
        raise ArithmeticError('No positive curvature after the uniform comparison loss.')
    base = finite.lower()+2*theta.lower()-arb(fcredit)-arb(constant_loss)
    gain = freeze_upper(arb(slope_cap)**2/curvature)
    family_floor = base-arb(gain)
    fixed_floor = arb(clin)+2*theta.lower()-arb(positive_far)-arb(old_error)-arb(fcredit)
    scalars = dict(reference_finite_expression=finite,
        reference_finite_derivative=derivative, negative_rho_energy_upper=arb(negative_rho_up),
        negative_far_square_credit_upper=arb(fcredit),
        positive_derivative_cap=phi_positive, negative_derivative_cap=phi_negative,
        signed_derivative_absolute_cap=arb(phi_cap), total_slope_absolute_cap=arb(slope_cap),
        lower_comparison_relative_error_upper=eps,
        uniform_comparison_constant_loss_upper=arb(constant_loss),
        uniform_comparison_curvature_loss_upper=arb(curvature_loss),
        comparison_curvature_lower=curvature, comparison_base_lower=base,
        maximal_steering_gain_upper=arb(gain), all_alpha_comparison_floor_lower=family_floor,
        fixed_envelope_selfterm_refinement_floor_lower=fixed_floor)
    require_finite(list(scalars.values()))
    exact = dict(alpha0=alpha, legal_coefficient=t0, frozen_positive_fee=fee,
        negative_rho_energy_upper=negative_rho_up, negative_far_square_credit_upper=fcredit,
        signed_derivative_absolute_cap=phi_cap, total_slope_absolute_cap=slope_cap,
        uniform_comparison_constant_loss_upper=constant_loss,
        uniform_comparison_curvature_loss_upper=curvature_loss,
        maximal_steering_gain_upper=gain, frozen_linear_envelope_upper=clin)
    return scalars, exact


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
            raise ValueError('Different accepted steering input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw))
    old, mixed, signed = [packets[name] for name in PINS]
    inherited = dict(PINS)
    for field in ('invocation_source_sha256', 'inherited_dependency_sha256'):
        for name, digest in signed[field].items():
            if name in inherited and inherited[name] != digest:
                raise ValueError('Conflicting saved steering dependency: '+name)
            inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked steering dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if ([p['schema'] for p in (old, mixed, signed)] != [
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1',
            'ACTUAL_C27_COHERENT_MIXED_CAPACITY_V1',
            'ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any(signed[k] != old[k] or old[k] != mixed[k] for k in
                   ('source_integer_sha256', 'source_bits', 'source_row'))
            or (signed['source_bits'], signed['source_row']) != (384, 0)
            or old['remaining_modes'] != mixed['remaining_modes']
            or len(old['remaining_modes']) != 18
            or mixed['inherited_capacity_packet_sha256'] !=
                PINS['full-gamma-negative-channel-capacity.json.zlib']):
        raise ValueError('Different same-source, core, lift, scalar grid or packet roles.')
    if any(mixed['exact_rationals'][k] != value for k, value in old['exact_rationals'].items()):
        raise ValueError('Different fixed legal coefficient, fee or steering origin.')
    rational = lambda packet, key: fmpq(*map(int, packet['exact_rationals'][key]))
    if rational(mixed, 'frozen_linear_envelope_upper') != rational(signed, 'frozen_linear_envelope_upper'):
        raise ValueError('Different unchanged linear envelope.')
    b = lambda packet, key: interval_ball(packet['scalar_bounds'][key], BITS)
    alpha, t0, fee = [rational(old, key) for key in
                     ('fixed_alpha', 'legal_coefficient', 'positive_source_valid_fee')]
    delta = arb(b(old, 'delta_lower').lower()).union(b(old, 'delta_upper').upper())
    scalars, exact = steering_floor(alpha, t0, fee, b(old, 'selected_reference_form'),
        b(old, 'finite_negative_channel'), delta, b(old, 'steering_H11_reference'),
        b(old, 'steering_H12_reference'), b(old, 'steering_reference_rRr').upper(),
        b(mixed, 'reference_far_energy_upper').upper(),
        b(old, 'negative_far_channel_radius_upper').upper(),
        b(signed, 'actual_signed_theta'), b(old, 'neumann_eta_upper').upper(),
        rational(mixed, 'frozen_linear_envelope_upper'),
        b(mixed, 'combined_reference_error_upper').upper())
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    if any(pair is None or int(pair[0]) > int(pair[1]) for pair in exported.values()):
        raise ArithmeticError('Invalid steering scalar export; stop without retry.')
    excluded = scalars['all_alpha_comparison_floor_lower'].lower() > 0
    result = dict(schema='ACTUAL_C27_ALL_ALPHA_COMPARISON_CAPACITY_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: digest for p, digest in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: digest for p, digest in guards.items() if p not in invoked},
        source_integer_sha256=old['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=1024, scalar_bounds=exported,
        exact_rationals={k: [str(v.p), str(v.q)] for k, v in exact.items()},
        comparison_domain='all real alpha, SAME fixed source/core/lift, frozen positive fee and legal t0',
        all_alpha_frozen_fee_comparison_excluded=bool(excluded),
        verdict='fixed-fee-all-alpha-comparison-excluded' if excluded else 'steering-capacity-not-excluded',
        actual_center_sign_claim=False, actual_full_border_sign_claim=False, actual_rh_proof=False,
        new_scalar_capacity_application=True, new_signed_derivative_suppliers=0,
        new_source_feature_matrix_inverse_trace_or_tail_calls=0, old_producers_replayed=0,
        source_parameter_precision_scan=False, elapsed_seconds=time.monotonic()-started,
        conditions=signed['conditions']+[
            'The affine row y(alpha)=alpha*rho+betaJ and the actual far row f use one fixed common source, old core, lift, fee and legal t0 for every real alpha.',
            'Qref=T R T*-e0e0*/delta, h11=1/t0+rho Qref rho*, h12=rho Qref betaJ*, r0=y0 T R T* y0* and theta0=f Qref y0* are the SAME reference objects.',
            'The saved inverse comparison R/(1+eta)<=S^-1 supplies Qtrue>=Qref-eta/(1+eta)*T R T*. All derivative and norm caps apply to that same realization.',
            'The floor concerns E_D(alpha)=alpha^2/t0+D_frozen+(y(alpha)+f)Qtrue(y(alpha)+f)*. A frozen fee upper cannot give a lower bound on the actual smaller-fee center.'],
        scope='A scalar-only capacity test for one fixed-fee/legal-coefficient all-alpha comparison family. Not a lower bound for the actual center, an RH refutation, a new generic theorem, Lean proof, full-border or cofinal result.')
    verify_guards(guards)
    digest = write_packet(out, result)
    print('saved', out.name, digest, flush=True)
    for key in ('negative_far_square_credit_upper', 'signed_derivative_absolute_cap',
                'comparison_curvature_lower', 'maximal_steering_gain_upper',
                'all_alpha_comparison_floor_lower'):
        print(key, scalars[key].str(24), flush=True)
    print(result['verdict'], 'seconds', time.monotonic()-started, flush=True)


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
