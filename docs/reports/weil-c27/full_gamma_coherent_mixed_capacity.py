"""Price one same-source coherent mixed direction before supplying its tail.

Reuse the frozen negative-channel packet, true core, exact17-reference,
full paid far Gram and source/trace bounds. One reference-vector action
and one18-dimensional Gram quadratic test a different linear envelope.
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
from full_gamma_negative_channel_capacity import PINS as BASE_PINS, SOURCE_ID
from joint_gamma_boundary import matrix_in, matrix_out
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
CAPACITY_NAME = 'full-gamma-negative-channel-capacity.json.zlib'
CAPACITY_SHA = '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e'
CAPACITY_SOURCE_SHA = 'f5ea417fa0aafdbff7e69a1d0305ce806b257290c439947a245976817995fdd1'
PINS = {**BASE_PINS, CAPACITY_NAME: CAPACITY_SHA}


def capacity(sigma, reference, rho, beta, gram, alpha, t0, fee, energy,
             kappa, eta, raw_r_upper, cplus_upper, lam, delta_lower, delta_upper):
    """Enclose one exact reference direction, not an independently chosen center."""
    size = sigma.nrows()
    if (size < 2 or sigma.ncols() != size
            or (reference.nrows(), reference.ncols()) != (size-1, size-1)
            or (gram.nrows(), gram.ncols()) != (size, size)
            or (rho.nrows(), rho.ncols(), beta.nrows(), beta.ncols()) != (1, size, 1, size)):
        raise ValueError('Incompatible same-core, reference, rows or full Gram.')
    values = (list(sigma.entries())+list(reference.entries())+list(gram.entries())
              +list(rho.entries())+list(beta.entries())
              +[arb(v) for v in (alpha, t0, fee, energy, kappa, eta, raw_r_upper,
                                 cplus_upper, lam, delta_lower, delta_upper)])
    if any(not v.is_finite() for v in values):
        raise ArithmeticError('Nonfinite input; stop without retry.')
    a, delta = sigma[0, 0], -sigma[0, 0]
    if (not a < 0 or not delta.lower() > 0 or not 0 < delta_lower <= delta_upper
            or not 0 <= eta < 1 or not t0 > 0 or not fee > 0 or not energy >= 0
            or not kappa > 0 or not raw_r_upper >= 0 or not cplus_upper >= 0):
        raise ArithmeticError('Missing accepted pivot, positive parameter or reference bound.')
    T = arb_mat([[-sigma[0, j+1]/a if i == 0 else arb(int(i == j+1))
                  for j in range(size-1)] for i in range(size)])
    y = arb(alpha)*rho+beta
    g = y*T
    action = reference*g.transpose()
    direction = T*action
    direction[0, 0] -= lam/delta
    if any(not v.is_finite() for v in direction.entries()):
        raise ArithmeticError('Nonfinite reference direction; stop without retry.')
    # r is the CACHED raw gRg*, never the already-inflated positive-energy bound.
    Aref = ((1+eta)*cplus_upper**2).upper()
    relative = (eta/(1-eta)).upper()
    error = (relative*(raw_r_upper.sqrt().upper()+Aref.sqrt().upper())**2).upper()
    lambda_lower = max(arb(0), abs(lam).lower())
    total = (arb(alpha)**2/arb(t0)+arb(fee)+raw_r_upper+Aref+error
             -lambda_lower**2/delta_upper)
    pair = envelope(total, BITS)
    if pair is None:
        raise ArithmeticError('Nonfinite linear-envelope upper; stop without retry.')
    frozen = fmpq(int(pair[1]), 1 << BITS)
    quadratic = (direction.transpose()*gram*direction)[0, 0]
    if not quadratic.is_finite() or not quadratic.upper() >= 0:
        raise ArithmeticError('Invalid full-Gram quadratic upper; stop without retry.')
    radius = (arb(kappa)*(arb(energy)*quadratic.upper()).sqrt()).upper()
    floor = (arb(frozen)-2*radius).lower()
    uniform_upper = (arb(frozen)+2*radius).upper()
    scalars = dict(neumann_eta_upper=eta, cached_raw_reference_energy_upper=raw_r_upper,
        cached_positive_far_norm_upper=cplus_upper, finite_negative_channel=lam,
        cached_delta_lower=delta_lower, cached_delta_upper=delta_upper,
        finite_negative_channel_magnitude_lower=lambda_lower,
        reference_far_energy_upper=Aref, reference_relative_error_upper=relative,
        combined_reference_error_upper=error, linear_envelope_before_export=total,
        frozen_linear_envelope_upper=arb(frozen), coherent_full_far_gram_quadratic=quadratic,
        coherent_full_far_gram_quadratic_upper=quadratic.upper(),
        coherent_response_radius_upper=radius, linear_package_floor_lower=floor,
        uniform_center_upper=uniform_upper, strict_theta_upper_threshold=arb(-frozen/2))
    if any(not v.is_finite() for v in scalars.values()):
        raise ArithmeticError('Nonfinite scalar output; stop without retry.')
    return direction, frozen, scalars


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            invoked.add(Path(filename).resolve())
    guards = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    needed = {'joint-gamma-6765-rho32-v2-schur.json', 'full-gamma-directional-metric.json.zlib',
              'full-gamma-lift-projection.json.zlib', 'full-gamma-projection-upper.json.zlib',
              'full-gamma-actual-channel-projection.json.zlib', CAPACITY_NAME}
    packets = {}
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted packet: '+name)
        if name in needed:
            packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    saved = packets[CAPACITY_NAME]
    if (saved['schema'] != 'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1'
            or saved['source_sha256'] != CAPACITY_SOURCE_SHA
            or saved['packet_sha256'] != BASE_PINS or saved['grid_bits'] != BITS
            or (saved['source_bits'], saved['source_row'], saved['source_integer_sha256'])
               != (384, 0, SOURCE_ID)):
        raise ValueError('Different frozen source, parameters or inherited certificate chain.')
    # Accepted identities, inertia, Neumann, trace and source tests are inherited.
    producers = {}
    inherited_sources = dict(saved['cached_producer_source_sha256'])
    inherited_sources['full_gamma_negative_channel_capacity.py'] = CAPACITY_SOURCE_SHA
    for name, digest in inherited_sources.items():
        path = directory/name
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Different cached producer source: '+name)
        producers[path] = digest
    core = packets['joint-gamma-6765-rho32-v2-schur.json']
    cert = packets['full-gamma-directional-metric.json.zlib']
    lift = packets['full-gamma-lift-projection.json.zlib']
    projection = packets['full-gamma-projection-upper.json.zlib']
    actual = packets['full-gamma-actual-channel-projection.json.zlib']
    modes = saved['remaining_modes']
    if len(modes) != 18 or any(p['remaining_modes'] != modes for p in
                             (core, cert, lift, projection, actual)):
        raise ValueError('Different same-core lift coordinate order.')
    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], BITS))
    reference = arb_mat([[arb(fmpq(v)) for v in row] for row in cert['reference_inverse']])
    rho = matrix_in([projection['matrices']['source_on_lift'][0]], BITS)
    beta = matrix_in(actual['actual_finite_row_on_lift'], BITS)
    gram = symmetric(matrix_in(lift['far_gram_upper'], BITS))
    exact = {k: fmpq(*map(int, saved['exact_rationals'][k])) for k in
             ('fixed_alpha', 'legal_coefficient', 'positive_source_valid_fee',
              'actual_source_energy_cap', 'kappa')}
    b = saved['scalar_bounds']
    upper = lambda key: interval_ball(b[key], BITS).upper()
    direction, frozen, scalars = capacity(sigma, reference, rho, beta, gram,
        exact['fixed_alpha'], exact['legal_coefficient'], exact['positive_source_valid_fee'],
        exact['actual_source_energy_cap'], exact['kappa'], upper('neumann_eta_upper'),
        upper('selected_reference_form'), upper('positive_far_channel_norm_upper'),
        interval_ball(b['finite_negative_channel'], BITS),
        interval_ball(b['delta_lower'], BITS).lower(), upper('delta_upper'))
    exports = {key: envelope(value, BITS) for key, value in scalars.items()}
    coefficients = matrix_out(direction, BITS)
    if (any(p is None or int(p[0]) > int(p[1]) for p in exports.values())
            or any(p is None or int(p[0]) > int(p[1]) for row in coefficients for p in row)):
        raise ArithmeticError('Invalid ordered outward export; stop without retry.')
    stopped = scalars['linear_package_floor_lower'] >= 0
    uniform_negative = scalars['uniform_center_upper'] < 0
    verdict = ('uniform-negative-center-upper' if uniform_negative else
               'fixed-coherent-linear-envelope-excluded' if stopped else
               'fixed-coherent-linear-envelope-not-excluded')
    result = dict(schema='ACTUAL_C27_COHERENT_MIXED_CAPACITY_V1', packet_sha256=PINS,
        inherited_capacity_packet_sha256=CAPACITY_SHA, input_sha256=saved['input_sha256'],
        source_sha256=guards[source], invocation_source_sha256={p.name: h for p, h in guards.items()},
        cached_producer_source_sha256={p.name: h for p, h in producers.items()},
        remaining_modes=modes, source_integer_sha256=SOURCE_ID, source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=1024, exact_rationals=saved['exact_rationals'] | {
            'frozen_linear_envelope_upper': [str(frozen.p), str(frozen.q)],
            'strict_theta_upper_threshold': [str((-frozen/2).p), str((-frozen/2).q)]},
        reference_direction_on_lift=coefficients, scalar_bounds=exports, verdict=verdict,
        actual_signed_theta_supplied=False, uniform_center_upper_negative=bool(uniform_negative),
        actual_full_border_sign_claim=False, steering_or_direction_optimality_claim=False,
        new_reference_vector_actions=1, new_reference_quadratic_forms=0,
        new_full_far_gram_forms=1, new_inverse_or_solve_calls=0, new_trace_contractions=0,
        new_source_or_far_suppliers=0, old_producers_replayed=0, new_determinant_calls=0,
        parameter_source_precision_scan=False, elapsed_seconds=time.monotonic()-started,
        conditions=saved['conditions']+[
            'The cached raw selected_reference_form is r=(yT)R(yT)* for SAME fixed alpha,rho,beta,T,R; the inflated positive energy is not substituted.',
            'vR=T*R*(yT)*-(lambda/delta)*e0 encloses ONE exact reference direction selected only from finite cached data; it is not a free midpoint or the actual inverse action.',
            'R<=(1+eta)*S^-1 and S^-1-R<=eta*R/(1-eta) price the full combined positive-complement reference error; no reference=actual inverse identity is used.',
            'The SAME full paid Gram including all allowance bounds theta=f*vR for f=kappa*u_far*W_far*R_far*V, with the original source energy upper.',
            'Clin uses the lambda absolute LOWER and delta UPPER; the interval direction uses full lambda/delta enclosures from the same true core.',
            'The negative tail self-square is discarded. Nonnegative DOWN(Clin-2*Ctheta) excludes only this fixed linear upper envelope; a negative floor supplies no sign or favorable realization.',
            'A SAME-source signed interval Itheta must satisfy UP(Clin+2*Itheta)<0; saved direction coefficient intervals must be propagated, not replaced by their midpoint for free.',
            'Only a negative actual center upper permits the additional PSD1263-good-block metric test; that metric is unpaid here.'],
        scope='Conditional cached-data application of existing block inverse, relative Loewner '
              'and PSD Cauchy bounds. No new generic method or Lean result. An actual signed '
              'tail interval, full-border sign, good-block metric, retained/cofinal positivity, '
              'Robin tail and RH remain unresolved.')
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
    print(verdict, 'seconds', result['elapsed_seconds'], flush=True)


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
