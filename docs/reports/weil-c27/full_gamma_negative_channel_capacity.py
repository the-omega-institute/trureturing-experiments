"""Test one fixed full-border envelope before supplying a signed far moment.

Reuse the TRUE old core, cached17-reference inverse and complete far trace.
Three small forms choose one tail-independent dyadic direction and price
its positive channel. The stopping floor excludes only this envelope.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, arb_mat, ctx, fmpq

from diagonal_floor_witness import interval_ball
from full_gamma_correlated_sign import symmetric
from joint_gamma_boundary import matrix_in
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS, ALPHA_BITS = 640, 64
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
PINS = {
    'joint-gamma-6765-rho32-v2-schur.json':
        '13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd',
    'full-gamma-directional-metric.json.zlib':
        'f6bea1b9b43c5ca8c85fc6cd7222038f2a40b81ed1ce2ef873af91ce1ca987a2',
    'full-gamma-reference-trace.json':
        'e6ea65765b341fef97051b6c484153526253452de78f135d8d570bfc411e57a1',
    'full-gamma-lift-projection.json.zlib':
        'a41813b27cd4567af326d42ec77e3b8c8942474b80bdf5833f7bc05aff52f873',
    'full-gamma-projection-upper.json.zlib':
        'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3',
    'full-gamma-actual-channel-projection.json.zlib':
        '63e3c2e291940ca3549fa4a1303b55b4d92c08d594a50daf7ca2132087ad2745',
    'full-gamma-common-column.json.zlib':
        '314e341189f72516aa9f251c6b873587f1bfae54d902c35f6c87c9b03130dfe7',
    'full-gamma-hankel-tail-energy.json.zlib':
        'd09a5ab254d164cfd40982b31e1a155e31d6288b705d916bf248236a0d527606',
    'full-gamma-surplus-budget.json':
        '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f',
    'full-gamma-correlated-sign.json':
        'aedb0a36cb4cec1f862ca0e41e64f76d39f10bd6da156e8404c10e04c4005149',
}
PRODUCERS = ('structured_residual_schur_finish.py', 'full_gamma_directional_metric.py',
    'full_gamma_reference_trace.py', 'full_gamma_lift_projection.py',
    'full_gamma_projection_upper.py', 'full_gamma_actual_channel_projection.py',
    'full_gamma_common_column.py', 'full_gamma_hankel_tail_energy.py',
    'full_gamma_surplus_budget.py', 'full_gamma_correlated_sign.py')


def capacity(sigma, reference, rho, beta, eta, t0, fee, energy,
             old_cap_lower, far_trace_upper, g00_upper):
    """Real old-core application; no inverse, solve, trace or far action."""
    size = sigma.nrows()
    if (size < 2 or sigma.ncols() != size
            or (reference.nrows(), reference.ncols()) != (size-1, size-1)
            or (rho.nrows(), rho.ncols(), beta.nrows(), beta.ncols()) != (1, size, 1, size)):
        raise ValueError('Incompatible true-core, reference or real row dimensions.')
    values = list(sigma.entries())+list(reference.entries())+list(rho.entries())+list(beta.entries())
    values += [eta, arb(t0), arb(fee), arb(energy), old_cap_lower, far_trace_upper, g00_upper]
    if any(not v.is_finite() for v in values):
        raise ArithmeticError('Nonfinite input; stop without retry.')
    a, delta = sigma[0, 0], -sigma[0, 0]
    if (not a < 0 or not delta.lower() > 0 or not 0 <= eta < 1
            or not t0 > 0 or not fee > 0 or not energy >= 0
            or not old_cap_lower > 0 or not far_trace_upper >= 0 or not g00_upper >= 0):
        raise ArithmeticError('Missing pivot, positivity or reference certificate guard.')
    T = arb_mat([[-sigma[0, j+1]/a if i == 0 else arb(int(i == j+1))
                  for j in range(size-1)] for i in range(size)])
    rT, bT = rho*T, beta*T
    rr = (rT*reference*rT.transpose())[0, 0]
    rb = (rT*reference*bT.transpose())[0, 0]
    h11 = 1/arb(t0)+rho[0, 0]**2/a+rr
    h12 = rho[0, 0]*beta[0, 0]/a+rb
    pair11, pair12 = envelope(h11, BITS), envelope(h12, BITS)
    if pair11 is None or pair12 is None or not h11.lower() > 0:
        raise ArithmeticError('No finite positive steering denominator; no fallback.')
    mid11, mid12 = (fmpq(int(p[0])+int(p[1]), 1 << (BITS+1)) for p in (pair11, pair12))
    ratio = -mid12/mid11
    alpha = fmpq((int(ratio.p)*(1 << ALPHA_BITS))//int(ratio.q), 1 << ALPHA_BITS)
    g = arb(alpha)*rT+bT
    gRg = (g*reference*g.transpose())[0, 0]
    if not gRg.is_finite() or not gRg.upper() >= 0:
        raise ArithmeticError('Invalid selected positive-form upper; stop without retry.')
    p_up = (gRg.upper()/(1-eta)).upper()
    kappa = arb(fmpq(33, 32))
    trace_up = (far_trace_upper/(kappa**2*old_cap_lower)).upper()
    cplus = (kappa*(arb(energy)*trace_up).sqrt()).upper()
    cminus = (kappa*(arb(energy)*g00_upper).sqrt()).upper()
    lam = arb(alpha)*rho[0, 0]+beta[0, 0]
    total = arb(alpha)**2/arb(t0)+arb(fee)+(p_up.sqrt()+cplus)**2
    total_pair = envelope(total, BITS)
    if total_pair is None:
        raise ArithmeticError('Nonfinite center-envelope upper; stop without retry.')
    frozen_upper = fmpq(int(total_pair[1]), 1 << BITS)
    floor = (arb(frozen_upper)-(abs(lam).upper()+cminus)**2/delta.lower()).lower()
    response_threshold = (delta.upper()*arb(frozen_upper)).sqrt().upper()
    scalars = dict(old_first_pivot=a, delta_lower=delta.lower(), delta_upper=delta.upper(),
        neumann_eta_upper=eta, steering_reference_rRr=rr, steering_reference_rRb=rb,
        steering_H11_reference=h11, steering_H12_reference=h12,
        selected_reference_form=gRg, positive_channel_energy_upper=p_up,
        old_trace_cap_constant_lower=old_cap_lower, reused_far_trace_upper=far_trace_upper,
        far_G00_upper=g00_upper, far_positive_channel_trace_upper=trace_up,
        positive_far_channel_norm_upper=cplus, negative_far_channel_radius_upper=cminus,
        finite_negative_channel=lam, center_envelope_before_export=total,
        frozen_center_envelope_upper=arb(frozen_upper), package_floor_lower=floor,
        sufficient_response_magnitude_upper=response_threshold)
    if any(not v.is_finite() for v in scalars.values()):
        raise ArithmeticError('Nonfinite scalar output; stop without retry.')
    return alpha, frozen_upper, scalars


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            invoked.add(Path(filename).resolve())
    guards = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    packets = []
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted packet: '+name)
        packets.append(json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw))
    core, cert, trace, lift, projection, actual, common, energy, budget, shared = packets
    if [p['schema'] for p in packets] != [
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1', 'ACTUAL_C27_FULL_GAMMA_DIRECTIONAL_METRIC_V1',
            'ACTUAL_C27_FULL_GAMMA_REFERENCE_TRACE_V1', 'ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
            'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1',
            'ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
            'ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1', 'ACTUAL_C27_SIGNED_HANKEL_TAIL_ENERGY_V1',
            'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1', 'ACTUAL_C27_FULL_GAMMA_CORRELATED_SIGN_V1']:
        raise ValueError('Different cached roles.')
    ch = PINS['joint-gamma-6765-rho32-v2-schur.json']
    sh = PINS['full-gamma-correlated-sign.json']
    if (cert['packet_sha256']['core'] != ch or trace['packet_sha256']['core'] != ch
            or shared['packet_sha256']['core'] != ch
            or trace['packet_sha256']['certificate'] != PINS['full-gamma-directional-metric.json.zlib']
            or cert['packet_sha256']['shared'] != sh or trace['packet_sha256']['shared'] != sh
            or cert['inherited_packet_sha256'] != shared['packet_sha256']
            or cert['verdict'] != 'directional-center-enclosed'
            or not cert['exact_reference_identities_checked'] or cert['positive_dimension'] != 17):
        raise ValueError('Different true core, accepted reference or trace scale.')
    for role, name in (('projection', 'full-gamma-lift-projection.json.zlib'),
                       ('common', 'full-gamma-common-column.json.zlib')):
        if shared['packet_sha256'][role] != PINS[name] or actual['packet_sha256'][role] != PINS[name]:
            raise ValueError('Different same-lift full far Gram or common source.')
    if actual['packet_sha256']['shared'] != sh:
        raise ValueError('Different finite row and full far trace correspondence.')
    for name in ('full-gamma-common-column.json.zlib',
                 'full-gamma-actual-channel-projection.json.zlib', 'full-gamma-surplus-budget.json'):
        if energy['packet_sha256'][name] != PINS[name]:
            raise ValueError('Different source-valid energy/fee/cap: '+name)
    if budget['packet_sha256']['full-gamma-common-column.json.zlib'] != PINS['full-gamma-common-column.json.zlib']:
        raise ValueError('Different h endpoint source.')
    modes, expected = core['remaining_modes'], core['input_sha256']
    if (len(modes) != 18 or modes[0] != 0 or any(
            p['remaining_modes'] != modes or p['input_sha256'] != expected
            for p in (cert, trace, lift, projection, actual, shared))
            or common['metadata'] != energy['metadata']
            or common['metadata']['input_sha256'] != expected
            or any(p['grid_bits'] != BITS for p in
                   (core, cert, trace, lift, projection, actual, common, energy, shared))):
        raise ValueError('Different basis, grid, lift order or scientific inputs.')
    if (projection['exact_source_bits'][0] != 384
            or projection['independent_source_indices'][0] != 0
            or projection['exact_source_integer_sha256'][0] != SOURCE_ID
            or any((p['source_bits'], p['source_row'], p['source_integer_sha256'])
                   != (384, 0, SOURCE_ID) for p in (energy, budget))
            or common['metadata']['basis'] !=
                '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'
            or not any('full allowance' in s for s in lift['conditions'])):
        raise ValueError('Different original unphased source or complete far allowance.')
    producers = {}
    for packet, name in zip(packets, PRODUCERS):
        path = directory/name
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != packet['source_sha256']:
            raise ValueError('Different cached producer: '+name)
        producers[path] = digest
    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], BITS))
    reference = arb_mat([[arb(fmpq(v)) for v in row] for row in cert['reference_inverse']])
    rho = matrix_in([projection['matrices']['source_on_lift'][0]], BITS)
    beta = matrix_in(actual['actual_finite_row_on_lift'], BITS)
    cap = fmpq(*map(int, energy['exact_rationals']['supplied_parameter_cap']))
    hlo, hhi = (fmpq(*map(int, budget['exact_rationals'][k])) for k in ('h_lower', 'h_upper'))
    if not 0 < cap < 1 or not 0 < hlo <= hhi:
        raise ArithmeticError('Require the saved positive eligibility cap and ordered h endpoints.')
    t0 = cap/hhi
    fee, ecap = (fmpq(int(energy['scalar_bounds'][k][1]), 1 << BITS)
                 for k in ('source_valid_fee_upper', 'far_energy_upper'))
    alpha, frozen, scalars = capacity(sigma, reference, rho, beta,
        interval_ball(cert['scalar_bounds']['neumann_eta_upper'], BITS).upper(), t0, fee, ecap,
        interval_ball(common['bounds']['far_action_energy_upper'], BITS).lower(),
        interval_ball(trace['scalar_bounds']['far_trace_upper'], BITS).upper(),
        interval_ball(lift['far_gram_upper'][0][0], BITS).upper())
    exported = {key: envelope(v, BITS) for key, v in scalars.items()}
    if any(p is None or int(p[0]) > int(p[1]) for p in exported.values()):
        raise ArithmeticError('Invalid ordered outward export; stop without retry.')
    stopped = scalars['package_floor_lower'] >= 0
    result = dict(schema='ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={p.name: h for p, h in guards.items()},
        cached_producer_source_sha256={p.name: h for p, h in producers.items()},
        input_sha256=expected, remaining_modes=modes, source_integer_sha256=SOURCE_ID,
        source_bits=384, source_row=0, grid_bits=BITS, working_bits=1024, alpha_grid_bits=ALPHA_BITS,
        exact_rationals={key: [str(v.p), str(v.q)] for key, v in dict(
            supplied_parameter_cap=cap, h_lower=hlo, h_upper=hhi, legal_coefficient=t0,
            positive_source_valid_fee=fee, actual_source_energy_cap=ecap,
            fixed_alpha=alpha, frozen_center_envelope_upper=frozen, kappa=fmpq(33, 32)).items()},
        scalar_bounds=exported, fixed_negative_channel_coordinate=0,
        verdict=('fixed-positive-channel-envelope-excluded' if stopped
                 else 'fixed-positive-channel-envelope-not-excluded'),
        actual_full_border_sign_claim=False, steering_optimality_claim=False,
        new_reference_quadratic_or_bilinear_forms=3, reference_form_dimension=17,
        new_trace_or_gram_contractions=0, new_inverse_or_solve_calls=0,
        new_source_or_far_suppliers=0, old_producers_replayed=0, new_determinant_calls=0,
        parameter_source_precision_scan=False, elapsed_seconds=time.monotonic()-started,
        conditions=energy['conditions']+trace['conditions']+[
            'The TRUE old core has its inherited one-negative inertia; Qplus=T*S^-1*T* and delta=-a use this SAME core, not SigmaMax.',
            'R/(1+eta)<=S^-1<=R/(1-eta) applies to the cached reference and actual positive17complement; no reference=actual identity is asserted.',
            'The original full far Gram, including the entire allowance, bounds the SAME R_far V; Ct=kappa^2*EM0*G uses the old common cap constant.',
            't0=cap/h_upper is a legal lower coefficient; D is the frozen positive source-valid fee. No shifted source row or tau*h fee is introduced.',
            'Alpha is one64-bit dyadic chosen by the finite reference steering ratio before any actual signed tail moment is supplied.',
            'A nonnegative package floor rules out only this fixed alpha and positive-channel envelope. A negative floor is inconclusive; favorable separate endpoints need not be jointly attained.',
            'A future SAME-source interval Ixi must give UP(Tbar-dist(0,Ilambda+Ixi)^2/delta_upper)<0 before the additional PSD good-block metric can be budgeted.'],
        scope='Conditional cached-data application of block inverse, certified reference Loewner, '
            'full far trace and Cauchy/Minkowski bounds. No new generic method or Lean result. '
            'Actual signed tail, other alpha/envelopes/parameters, full-border sign, good-block metric, '
            'retained/cofinal positivity, Robin tail and RH remain unresolved.')
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Accepted packet changed: '+name)
    for path, digest in {**guards, **producers}.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Invoked helper or cached producer changed: '+path.name)
    with out.open('xb') as stream:
        stream.write(zlib.compress(json.dumps(result, separators=(',', ':')).encode(), 9))
    for key, value in scalars.items():
        print(key, value.str(20), flush=True)
    print('fixed alpha', arb(alpha).str(20), 'legal t0', arb(t0).str(20), flush=True)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Use a fresh .zlib output; accepted archives are never overwritten.')
    ctx.prec = 1024
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
