"""Bound a same-source plane residual from an accepted raw-coordinate cap.

Reuse timelike projection and PSD Cauchy. Only saved scalars are combined;
no old scientific supplier or plane-floor calculation is repeated.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import zlib

from flint import arb, ctx, fmpq
from diagonal_floor_witness import interval_ball
from full_gamma_coherent_split_capacity import freeze_upper
from full_gamma_negative_index_joint import freeze_lower
from full_gamma_signed_mu import require_finite, verify_guards, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
PINS = {
    'full-gamma-cached-source-response.json.zlib':
        '2f6fa2f46f91dfa0449d5bec2523cdbd91a0e7587284f8e417ccb5bd2068ee03',
    'full-gamma-joint-plane-paid-floor.json.zlib':
        '3bdf5f3f81daeeda99c5fa71f2a48ee64c84c09da0f20070b45b0ad781aeec22',
    'full-gamma-relaxation-row-witness.json.zlib':
        '1879f4b0653e762ca0d7bb1e1933a59afc79fb10851d155f1c16c87055872992',
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
}


def coordinate_complement(b, g, c, t, a, gamma, lam, delta, theta, zeta, cap):
    """Inputs enclose ONE real form/row; lam and cap use raw coordinates."""
    require_finite([arb(x) for x in
                    (b, g, c, t, a, gamma, lam, delta, theta, zeta, cap)])
    if (not b > 0 or not a > 0 or not t >= 0 or not delta > 0
            or not 0 < gamma <= 1 or lam.contains(0) or cap < 0):
        raise ValueError('Missing common-plane, raw-coordinate or positive-form domain.')
    omega = c+g/b
    if not a.overlaps(t+b*omega*omega):
        raise ArithmeticError('Inconsistent common-plane norm; overlap proves no identity.')
    # xi=-delta*n, n in ker(P), n[0]=1, so Q*xi*=e0 in the SAME form.
    norm = -delta+lam*lam*t/(b*a)
    offset = lam*(t*theta/(b*a)-arb(gamma)*omega*zeta/a)
    lo = freeze_lower(offset.lower())-cap
    hi = freeze_upper(offset.upper())+cap
    nlo, nhi = freeze_lower(norm.lower()), freeze_upper(norm.upper())
    if lo > hi or nhi < 0:
        raise ArithmeticError('Inconsistent PSD complement or coordinate interval.')
    distance = max(fmpq(0), lo, -hi)
    if nhi == 0 and distance > 0:
        raise ArithmeticError('Nonzero cross on a null PSD-complement direction.')
    certified = nlo > 0
    lower = distance*distance/nhi if certified else fmpq(0)
    return dict(omega=omega, projected_coordinate_norm=norm,
                plane_coordinate_offset=offset), dict(
        cross_lower=lo, cross_upper=hi, cross_distance_lower=distance,
        projected_norm_lower=nlo, projected_norm_upper=nhi, residual_lower=lower), dict(
        positive_complement_norm_certified=bool(certified),
        coordinate_cross_separated=bool(distance > 0),
        positive_residual_lower_certified=bool(lower > 0))


def consume(out):
    started = time.perf_counter()
    source = Path(__file__).resolve()
    guards = {source: hashlib.sha256(source.read_bytes()).hexdigest()}
    packets = {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted coordinate-complement input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw))
    response, plane, relaxed, old = packets.values()
    inherited = dict(PINS)
    for packet in (response, plane, relaxed):
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256'):
            for name, digest in packet[field].items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting inherited input: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different current dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_CACHED_SOURCE_RESPONSE_V1', 'ACTUAL_C27_JOINT_PLANE_PAID_FLOOR_V1',
            'C27_SAVED_CONSTRAINT_RELAXATION_ROW_V1',
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any((p['source_bits'], p['source_row'], p['source_integer_sha256']) !=
                   (384, 0, response['source_integer_sha256']) for p in (response, plane, old))
            or (relaxed['reference_source_bits'], relaxed['reference_source_row'],
                relaxed['reference_source_integer_sha256']) !=
                (384, 0, response['source_integer_sha256'])
            or response['packet_sha256']['full-gamma-joint-plane-paid-floor.json.zlib'] !=
                PINS['full-gamma-joint-plane-paid-floor.json.zlib']
            or any(plane['packet_sha256'][n] != PINS[n] for n in tuple(PINS)[2:])
            or any(response[k] for k in
                   ('actual_center_sign_claim', 'actual_full_border_sign_claim', 'actual_rh_proof'))):
        raise ValueError('Different common source/form/coordinate contract.')
    rational = lambda p, key: fmpq(*map(int, p['exact_rationals'][key]))
    ball = lambda p, key: interval_ball(p['scalar_bounds'][key], BITS)
    if (rational(relaxed, 'reference_alpha0') != rational(old, 'fixed_alpha')
            or rational(plane, 'alpha0') != rational(old, 'fixed_alpha')
            or any(response['exact_rationals']['fixed_theta_'+end] !=
                   response['exact_rationals']['energy_endpoint_theta_'+end]
                   for end in ('lower', 'upper'))):
        raise ValueError('Different reference frame or transported response intervals.')
    # First actual new scalar application; no old floor is recalculated.
    zlo, zhi = [rational(response, 'zeta_'+end) for end in ('lower', 'upper')]
    if zlo > zhi:
        raise ValueError('Reversed accepted source-response interval.')
    zeta = arb(arb((zlo+zhi)/2), arb((zhi-zlo)/2))
    scalars, exact, flags = coordinate_complement(
        ball(plane, 'timelike_energy'), ball(plane, 'rho_timelike_pairing'),
        ball(plane, 'coordinate_ratio'), ball(plane, 'zero_coordinate_energy'),
        ball(response, 'fixed_positive_plane_energy'), rational(plane, 'gamma'),
        ball(old, 'finite_negative_channel'), -ball(old, 'old_first_pivot'),
        ball(response, 'fixed_transported_theta'), zeta,
        rational(relaxed, 'relaxed_negative_coordinate_cap'))
    recipes = {}
    for recipe in ('fixed', 'energy_endpoint'):
        same_f = rational(response, recipe+'_all_alpha_floor_lower')
        previous = rational(response, recipe+'_retained_floor_lower')
        augmented = same_f+exact['residual_lower']
        retained = max(previous, augmented)
        exact.update({recipe+'_same_F_floor_lower': same_f,
                      recipe+'_previous_retained_lower': previous,
                      recipe+'_residual_augmented_lower': augmented,
                      recipe+'_retained_lower': retained})
        recipes[recipe] = dict(strictly_improved=bool(retained > previous),
                               positive_paid_lower_certified=bool(retained > 0))
    exported = {key: envelope(v, BITS) for key, v in scalars.items()}
    if any(p is None or int(p[0]) > int(p[1]) for p in exported.values()):
        raise ArithmeticError('Invalid outward scalar export; no retry.')
    result = dict(schema='ACTUAL_C27_COORDINATE_COMPLEMENT_V1', packet_sha256=PINS,
        source_sha256=guards[source], inherited_dependency_sha256={
            p.name: h for p, h in guards.items() if p != source},
        source_integer_sha256=response['source_integer_sha256'],
        source_bits=384, source_row=0, grid_bits=BITS, working_bits=1024,
        scalar_bounds=exported,
        exact_rationals={key: [str(v.p), str(v.q)] for key, v in exact.items()},
        **flags, recipes=recipes, new_coordinate_complement_applications=1,
        new_paid_lower_additions=2, new_plane_floor_calls=0,
        new_source_matrix_inverse_solve_or_tail_suppliers=0, old_producers_or_fixtures_replayed=0,
        parameter_source_precision_cutoff_scan=False, actual_center_sign_claim=False,
        actual_full_border_sign_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started,
        conditions=[
            'All inherited actual-form, common-source, source-energy and physical-error contracts remain. P=T R T*, T has identity tail and R is positive definite; n=(1,sigma[0,1:]/sigma[0,0]) lies in ker(P).',
            'Q=gamma*P-e0e0*/delta, xi=-delta*n, Q*xi*=e0. Its projection off y0 and v lies in the SAME nonnegative complement. Lambda and cap use raw finite-row units, not division by sqrt(delta).',
            'The cap bounds the SAME actual |f[0]|. Offset and norm identities require shared parameters, not favorable independent endpoints or interval overlap as evidence of realization.',
            'Only certified positive projected norm permits Cauchy division. A zero-containing cross gives residual lower0. A norm interval not certified positive gives no positive gain here.',
            'Add this uniform residual lower only to the saved lower of the SAME F. Take a maximum with previous retained lowers afterward; no residual is added to an unrelated older bound.',
            'A nonpositive lower is inconclusive. No actual energy upper, negative physical row, actual-center/full-border/good-block/cofinal positivity, Robin or RH proof is supplied. Classical projection/Cauchy only; no novelty or Lean claim.',
        ])
    verify_guards(guards)
    if Path(out).exists():
        raise FileExistsError(out)
    digest = write_packet(out, result)
    print('saved coordinate-complement result', digest, 'seconds', result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    ctx.prec = 1024
    consume(args.out)
