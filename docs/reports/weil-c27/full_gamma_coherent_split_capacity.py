"""Bound the saved coherent direction by the already certified PSD trace.

Scalar-only consumption of two accepted packets. No matrix, vector, trace,
source or tail producer is evaluated. This tests one fixed upper envelope.
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
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
PINS = {
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
    'full-gamma-coherent-mixed-capacity.json.zlib':
        'b082eaad102c5f7ca38d1c49c19970952288ce630198a4ea76b48b96d8d23b65',
}
SOURCE_PINS = {
    'full_gamma_negative_channel_capacity.py':
        'f5ea417fa0aafdbff7e69a1d0305ce806b257290c439947a245976817995fdd1',
    'full_gamma_coherent_mixed_capacity.py':
        '81acfc63c4328b7e69a1827955a2985e880e5d03a626d890458e36a4aee5d166',
}


def freeze_upper(value):
    pair = envelope(value, BITS)
    if pair is None or int(pair[0]) > int(pair[1]):
        raise ArithmeticError('Nonfinite or unordered scalar export; stop without retry.')
    return fmpq(int(pair[1]), 1 << BITS)


def split_capacity(r, eta, trace, lam, delta_lower, g00, direct_q,
                   energy, kappa, clin):
    """Bound the SAME exact majorant seminorm, not a coefficient-box norm."""
    values = [arb(v) for v in (r, eta, trace, lam, delta_lower, g00,
                               direct_q, energy, kappa, clin)]
    if any(not v.is_finite() for v in values):
        raise ArithmeticError('Nonfinite scalar input; stop without retry.')
    if (not r >= 0 or not 0 <= eta < 1 or not trace >= 0
            or not delta_lower > 0 or not g00 >= 0 or not direct_q >= 0
            or not energy >= 0 or not kappa > 0):
        raise ArithmeticError('Missing certified scalar domain; stop without retry.')
    positive = (r*(1+eta)*trace).sqrt().upper()
    negative = (abs(lam).upper()/delta_lower*g00.sqrt().upper()).upper()
    split = freeze_upper(positive+negative)
    direct = freeze_upper(direct_q.sqrt().upper())
    cap = min(split, direct)
    radius = freeze_upper(arb(kappa)*arb(energy).sqrt().upper()*arb(cap))
    floor = (arb(clin)-2*arb(radius)).lower()
    uniform = (arb(clin)+2*arb(radius)).upper()
    scalars = dict(cached_raw_reference_energy_upper=r, cached_neumann_eta_upper=eta,
        cached_positive_channel_trace_upper=trace, finite_negative_channel=lam,
        cached_delta_lower=delta_lower, cached_far_G00_upper=g00,
        cached_direct_gram_quadratic_upper=direct_q,
        positive_direction_norm_upper=positive, negative_direction_norm_upper=negative,
        split_direction_norm_upper=arb(split), direct_direction_norm_upper=arb(direct),
        selected_direction_norm_upper=arb(cap), coherent_response_radius_upper=arb(radius),
        frozen_linear_envelope_upper=arb(clin), linear_package_floor_lower=floor,
        uniform_center_upper=uniform, strict_theta_upper_threshold=arb(-clin/2))
    if any(not v.is_finite() for v in scalars.values()):
        raise ArithmeticError('Nonfinite scalar result; stop without retry.')
    return scalars, cap, radius


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            invoked.add(Path(filename).resolve())
    guards = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    packets = {}
    for name, digest in PINS.items():
        path = directory/name
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted scalar packet: '+name)
        guards[path] = digest
        packets[name] = json.loads(zlib.decompress(raw))
    old = packets['full-gamma-negative-channel-capacity.json.zlib']
    saved = packets['full-gamma-coherent-mixed-capacity.json.zlib']
    if (old['schema'] != 'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1'
            or saved['schema'] != 'ACTUAL_C27_COHERENT_MIXED_CAPACITY_V1'
            or old['source_sha256'] != SOURCE_PINS['full_gamma_negative_channel_capacity.py']
            or saved['source_sha256'] != SOURCE_PINS['full_gamma_coherent_mixed_capacity.py']
            or saved['packet_sha256'] != old['packet_sha256'] | {
                'full-gamma-negative-channel-capacity.json.zlib':
                    PINS['full-gamma-negative-channel-capacity.json.zlib']}
            or saved['inherited_capacity_packet_sha256'] !=
                PINS['full-gamma-negative-channel-capacity.json.zlib']):
        raise ValueError('Different accepted source or inherited certificate chain.')
    for key in ('grid_bits', 'remaining_modes', 'input_sha256',
                'source_integer_sha256', 'source_bits', 'source_row'):
        if saved[key] != old[key]:
            raise ValueError('Different common source/lift metadata: '+key)
    if (saved['grid_bits'] != BITS or saved['source_bits'] != 384
            or saved['source_row'] != 0 or len(saved['remaining_modes']) != 18
            or any(saved['exact_rationals'][k] != v
                   for k, v in old['exact_rationals'].items())):
        raise ValueError('Different frozen fixed parameters, source or grid.')
    # Reuse certificate identities; only raw bytes are guarded, not replayed.
    inherited_sources = {}
    for packet in (old, saved):
        for field in ('cached_producer_source_sha256', 'invocation_source_sha256'):
            for name, digest in packet[field].items():
                if name in inherited_sources and inherited_sources[name] != digest:
                    raise ValueError('Conflicting inherited source pin: '+name)
                inherited_sources[name] = digest
    for name, digest in (saved['packet_sha256'] | inherited_sources | SOURCE_PINS).items():
        path = directory/name
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Different accepted dependency: '+name)
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked dependency: '+name)
        guards[path] = digest
    exact = {k: fmpq(*map(int, saved['exact_rationals'][k])) for k in
             ('frozen_linear_envelope_upper', 'strict_theta_upper_threshold',
              'actual_source_energy_cap', 'kappa')}
    clin = exact['frozen_linear_envelope_upper']
    if exact['strict_theta_upper_threshold'] != -clin/2:
        raise ValueError('Different frozen signed-response threshold.')
    def ball(packet, key):
        pair = packet['scalar_bounds'][key]
        if len(pair) != 2 or int(pair[0]) > int(pair[1]):
            raise ValueError('Unordered scalar input: '+key)
        return interval_ball(pair, BITS)
    upper = lambda packet, key: ball(packet, key).upper()
    scalars, norm, radius = split_capacity(
        upper(saved, 'cached_raw_reference_energy_upper'),
        upper(saved, 'neumann_eta_upper'),
        upper(old, 'far_positive_channel_trace_upper'),
        ball(saved, 'finite_negative_channel'),
        ball(saved, 'cached_delta_lower').lower(),
        upper(old, 'far_G00_upper'),
        upper(saved, 'coherent_full_far_gram_quadratic_upper'),
        exact['actual_source_energy_cap'], exact['kappa'], clin)
    exports = {k: envelope(v, BITS) for k, v in scalars.items()}
    if any(p is None or int(p[0]) > int(p[1]) for p in exports.values()):
        raise ArithmeticError('Nonfinite or unordered output; stop without retry.')
    stopped = scalars['linear_package_floor_lower'] >= 0
    uniform_negative = scalars['uniform_center_upper'] < 0
    verdict = ('uniform-negative-center-upper' if uniform_negative else
               'fixed-coherent-linear-envelope-excluded' if stopped else
               'fixed-coherent-linear-envelope-not-excluded')
    result = dict(schema='ACTUAL_C27_COHERENT_SPLIT_CAPACITY_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: h for p, h in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: h for p, h in guards.items() if p not in invoked},
        input_sha256=saved['input_sha256'], remaining_modes=saved['remaining_modes'],
        source_integer_sha256=saved['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=1024,
        exact_rationals=saved['exact_rationals'] | {
            'selected_direction_norm_upper': [str(norm.p), str(norm.q)],
            'coherent_response_radius_upper': [str(radius.p), str(radius.q)]},
        scalar_bounds=exports, verdict=verdict, actual_signed_theta_supplied=False,
        uniform_center_upper_negative=bool(uniform_negative), actual_full_border_sign_claim=False,
        new_reference_vector_actions=0, new_reference_quadratic_forms=0,
        new_full_far_gram_forms=0, new_trace_contractions=0, new_inverse_or_solve_calls=0,
        new_source_or_far_suppliers=0, old_producers_replayed=0, new_determinant_calls=0,
        parameter_source_precision_scan=False, elapsed_seconds=time.monotonic()-started,
        conditions=saved['conditions']+[
            'The SAME paid PSD G and R<=(1+eta)*S^-1 give tr(T R T* G)<=(1+eta)*Tfar. Tfar is already normalized; no energy division is repeated.',
            'PSD operator<=trace bounds the ONE analytic direction: Nsplit=sqrt(r*(1+eta)*Tfar)+abs(lambda)_UP/delta_LOWER*sqrt(G00).',
            'Ncap=min(Nsplit,Ndirect) bounds the SAME majorant seminorm, not all independently selected vectors in the saved coefficient box.',
            'The source energy is an upper cap. No particular-source functional bound is divided to infer a direction norm.',
            'A negative package floor supplies neither an actual signed response nor its favorable realization. No expensive supplier is warranted by nonexclusion alone.',
            'The exact Clin and signed threshold are unchanged. A future SAME-source theta interval must retain direction/lift relations and pay all errors before the extra PSD good-block metric.'
        ],
        scope='Scalar-only cached application of existing PSD trace, Loewner and Cauchy bounds; no generic method or Lean result. Actual theta, full-border sign, good-block metric, retained/cofinal positivity, Robin tail and RH remain unresolved.')
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise RuntimeError('Read dependency changed during scalar consumption: '+str(path))
    raw = json.dumps(result, separators=(',', ':'), allow_nan=False).encode()
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('xb') as stream:
        stream.write(zlib.compress(raw, 9))
    print(json.dumps({k: str(v) for k, v in scalars.items()}, indent=2))
    print(verdict)


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
