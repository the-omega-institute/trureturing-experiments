"""Upper-bound the SAME residual before paying for another source query.

Only accepted saved scalars are combined. This bounds the residual and
strength of the specified SAME-F lower algorithm, not actual paid energy.
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
    'full-gamma-joint-far-row.json.zlib':
        '00b9e4bdd8a99da4d6fdafc696d313bd43e4060e1ba7d8cec316f7cc5bcb7707',
    'full-gamma-cached-source-response.json.zlib':
        '2f6fa2f46f91dfa0449d5bec2523cdbd91a0e7587284f8e417ccb5bd2068ee03',
    'full-gamma-joint-plane-paid-floor.json.zlib':
        '3bdf5f3f81daeeda99c5fa71f2a48ee64c84c09da0f20070b45b0ad781aeec22',
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
}


def residual_upper_bound(b, g, c, a, gamma, theta, zeta, coordinate, delta,
                         positive_energy_upper):
    """Require the inherited common-form identities, not endpoint attainability."""
    require_finite([arb(v) for v in (b, g, c, a, gamma, theta, zeta,
                                    coordinate, delta, positive_energy_upper)])
    if (not b > 0 or not a > 0 or not delta > 0 or not 0 < gamma <= 1
            or not positive_energy_upper >= 0):
        raise ValueError('Missing common residual and positive-form domains.')
    if not isinstance(gamma, fmpq) or not isinstance(positive_energy_upper, fmpq):
        raise ValueError("Exact rational scale and positive-energy upper required.")
    omega = c+g/b
    response = arb(gamma)*zeta+omega*theta
    coordinate_distance = max(fmpq(0), freeze_lower(abs(coordinate).lower()))
    response_distance = max(fmpq(0), freeze_lower(abs(response).lower()))
    theta_cap = freeze_upper(abs(theta).upper())
    blo = freeze_lower(b.lower())
    aup, dup = freeze_upper(a.upper()), freeze_upper(delta.upper())
    if not blo > 0:
        raise ValueError('No positive directed timelike denominator.')
    upper = gamma*positive_energy_upper-coordinate_distance**2/dup \
        +theta_cap**2/blo-response_distance**2/aup
    if upper < 0:
        raise ArithmeticError('Residual upper contradicts its inherited nonnegativity.')
    exact = dict(residual_upper=upper, coordinate_distance_lower=coordinate_distance,
        plane_response_distance_lower=response_distance, theta_absolute_upper=theta_cap,
        timelike_energy_lower=blo, positive_plane_energy_upper=aup,
        negative_pivot_magnitude_upper=dup, far_positive_energy_upper=positive_energy_upper,
        gamma=gamma)
    return dict(plane_omega=omega, positive_plane_response=response), exact


def certificate_ceiling(same_F_lower, previous_lower, residual_upper):
    """Ceiling of an unchanged lower algorithm; no actual paid upper follows."""
    require_finite([arb(x) for x in (same_F_lower, previous_lower, residual_upper)])
    if any(not isinstance(x, fmpq) for x in (same_F_lower, previous_lower, residual_upper)):
        raise ValueError("Exact rational certificate inputs required.")
    if residual_upper < 0:
        raise ValueError('A nonnegative residual upper is required.')
    return max(previous_lower, same_F_lower+residual_upper)


def consume(out):
    started, source = time.perf_counter(), Path(__file__).resolve()
    guards = {source: hashlib.sha256(source.read_bytes()).hexdigest()}
    packets = {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted residual-envelope input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw))
    row, response, plane, old = packets.values()
    inherited = dict(PINS)
    for packet in packets.values():
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256',
                      'cached_producer_source_sha256'):
            for name, digest in packet.get(field, {}).items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting residual-envelope binding: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        p = HERE/name
        if p in guards and guards[p] != digest:
            raise ValueError('Different residual-envelope dependency: '+name)
        guards[p] = digest
    verify_guards(guards)
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_JOINT_FAR_ROW_V1', 'ACTUAL_C27_CACHED_SOURCE_RESPONSE_V1',
            'ACTUAL_C27_JOINT_PLANE_PAID_FLOOR_V1',
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any((p['source_bits'], p['source_row'], p['source_integer_sha256']) !=
                   (384, 0, row['source_integer_sha256']) for p in (response, plane, old))
            or any(row['packet_sha256'][n] != PINS[n] for n in tuple(PINS)[1:])
            or response['packet_sha256']['full-gamma-joint-plane-paid-floor.json.zlib'] !=
                PINS['full-gamma-joint-plane-paid-floor.json.zlib']
            or not row['reused_coordinate0'] or len(row['actual_far_row']) != 1
            or len(row['actual_far_row'][0]) != 18):
        raise ValueError('Different common source, raw coordinate or finite form.')
    rational = lambda p, key: fmpq(*map(int, p['exact_rationals'][key]))
    ball = lambda p, key: interval_ball(p['scalar_bounds'][key], BITS)
    if (any(response['exact_rationals']['fixed_theta_'+end] !=
            response['exact_rationals']['energy_endpoint_theta_'+end]
            for end in ('lower', 'upper'))
            or rational(plane, 'alpha0') != rational(old, 'fixed_alpha')):
        raise ValueError('Different SAME response or steering reference.')
    zlo, zhi = [rational(response, 'zeta_'+end) for end in ('lower', 'upper')]
    if zlo > zhi:
        raise ValueError('Reversed accepted response interval.')
    zeta = arb(arb((zlo+zhi)/2), arb((zhi-zlo)/2))
    scalars, exact = residual_upper_bound(
        ball(plane, 'timelike_energy'), ball(plane, 'rho_timelike_pairing'),
        ball(plane, 'coordinate_ratio'), ball(response, 'fixed_positive_plane_energy'),
        rational(plane, 'gamma'), ball(response, 'fixed_transported_theta'), zeta,
        interval_ball(row['actual_far_row'][0][0], BITS), -ball(old, 'old_first_pivot'),
        rational(plane, 'fixed_far_positive_energy_upper'))
    lower = rational(row, 'full_residual_lower')
    upper = min(exact['residual_upper'], rational(row, 'full_residual_upper'))
    if not 0 <= lower <= upper:
        raise ArithmeticError('New residual envelope conflicts with the saved lower.')
    exact.update(previous_residual_upper=rational(row, 'full_residual_upper'),
                 selected_residual_lower=lower, selected_residual_upper=upper)
    recipes = {}
    for recipe in ('fixed', 'energy_endpoint'):
        same_F = rational(response, recipe+'_all_alpha_floor_lower')
        previous = rational(row, recipe+'_retained_lower')
        ceiling = certificate_ceiling(same_F, previous, upper)
        exact.update({recipe+'_same_F_lower':same_F, recipe+'_previous_retained_lower':previous,
                      recipe+'_residual_only_retained_certificate_ceiling':ceiling})
        recipes[recipe] = dict(residual_only_refinement_cannot_certify_positive=bool(ceiling < 0))
    exported = {k:envelope(v, BITS) for k,v in scalars.items()}
    if any(v is None or int(v[0]) > int(v[1]) for v in exported.values()):
        raise ArithmeticError('Invalid residual-envelope export; no retry.')
    result = dict(schema='ACTUAL_C27_RESIDUAL_ENVELOPE_V1', packet_sha256=PINS,
        source_sha256=guards[source], inherited_dependency_sha256={p.name:h for p,h in guards.items() if p != source},
        source_integer_sha256=row['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=ctx.prec, scalar_bounds=exported,
        exact_rationals={k:[str(v.p),str(v.q)] for k,v in exact.items()}, recipes=recipes,
        new_residual_envelope_applications=1, new_certificate_ceilings=2,
        new_source_feature_matrix_reference_inverse_solve_or_tail_calls=0,
        old_scientific_producers_helpers_or_fixtures_replayed=0,
        source_direction_parameter_precision_cutoff_scan=False,
        actual_paid_energy_upper_claim=False, actual_center_sign_claim=False,
        actual_full_border_sign_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started, conditions=[
            'The SAME real Qeps=gamma*P-e0e0*/delta has one negative direction; y0 is timelike, v=rho+(g/b)y0 has positive squared norm a, and omega=c+g/b. All inherited common-source/form/physical-error premises remain.',
            'The actual fixed18-row has P(f,f)<=Aup and residual R=gamma*P(f,f)-f0^2/delta+theta^2/b-s^2/a>=0, s=gamma*zeta+omega*theta. Directed endpoint inequalities upper-bound this SAME residual; independent endpoint attainment is not asserted.',
            'The selected residual upper is the minimum of two upper bounds for the SAME residual, and must not contradict the previously positive lower. No lower certificates are recomputed or improved.',
            'With SAME saved F lower and previous retained lower fixed, max(previous,LF+Rup) bounds possible strength of residual-only LOWER refinement. It is not an upper for actual paid energy and does not exclude another F bound or a changed recipe.',
            'A nonnegative certificate ceiling is inconclusive. No new physical response, actual paid/center/full-border/cofinal positivity, Robin or RH proof. Standard projection and interval methods; no generic novelty or Lean claim.',
        ])
    verify_guards(guards)
    digest = write_packet(out,result)
    print('saved residual envelope',digest,'seconds',result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Require a fresh .zlib output; never overwrite accepted data.')
    ctx.prec = 1024
    consume(args.out)
