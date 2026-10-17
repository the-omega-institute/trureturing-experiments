"""Refine the saved zero-channel norm without replaying its finite actions.

Reuse the saved new-direction quadratic/features and existing directional
moment bound. The original sharp forward allowance is homogeneous in the
lift norm. No old or new direction, feature, quadratic or cross is repeated.
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
from full_gamma_directional_jet_norm import directional_norm
from full_gamma_negative_index_joint import freeze_lower
from full_gamma_shared_projection_energy import projected_energy
from full_gamma_signed_mu import require_finite, verify_guards, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS, R, M = 640, 56, 32768
PINS = {
    'full-gamma-zero-channel-pairing.json.zlib':
        'c4cbbead777f652812c0aa0bd238bf08053daa52be7afbedfa1321da19847a9b',
    'far-coupling.json.zlib':
        'ac35edacd7f55de470586e3b8a5915bcaa49b9d2c4c9c2ee073eeaaa8c9c73e7',
    'full-gamma-primal-tail.json':
        'fad634fb25d8d20e83a15d8277569978646246c2a560d0dc0f7a854b5087455a',
    'full-gamma-signed-correction.json.zlib':
        'd527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd',
}


def scaled_forward(old_forward, old_lift, new_lift):
    """Reuse the original uniform coefficient times its lift-norm cap."""
    require_finite([arb(x) for x in (old_forward, old_lift, new_lift)])
    if not old_lift > 0 or not old_forward >= 0 or not new_lift >= 0:
        raise ValueError('Missing homogeneous sharp-forward bound domains.')
    return freeze_upper(arb(old_forward)*arb(new_lift)/arb(old_lift))


def consume(out):
    started, source = time.monotonic(), Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name and Path(name).resolve().parent == HERE:
            invoked.add(Path(name).resolve())
    guards = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    packets = {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted zero-channel refinement input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    coarse, far, primal, signed = [packets[name] for name in PINS]
    inherited = dict(PINS)
    for packet in (coarse, signed):
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256'):
            for name, digest in packet[field].items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting saved refinement dependency: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked refinement dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_ZERO_CHANNEL_PAIRING_BOUND_V1', 'ACTUAL_C27_FAR_COUPLING_GRAM_V1',
            'ACTUAL_C27_FULL_GAMMA_PRIMAL_TAIL_V1', 'ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or coarse['source_integer_sha256'] != signed['source_integer_sha256']
            or (coarse['source_bits'], coarse['source_row']) != (384, 0)
            or (signed['source_bits'], signed['source_row']) != (384, 0)
            or (far['c'], far['N'], far['M'], far['jet_order'],
                far['moment_diagonals'], far['bound_grid_bits']) != (27, 1280, M, R, R, 628)
            or primal['metadata']['input_sha256'] != coarse['input_sha256']
            or primal['metadata']['basis'] !=
                '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'):
        raise ValueError('Different source, saved direction, physical moment mask or grids.')
    rational = lambda p, key: fmpq(*map(int, p['exact_rationals'][key]))
    ball = lambda p, key: interval_ball(p['scalar_bounds'][key], BITS)
    upper = lambda p, key: fmpq(int(p['scalar_bounds'][key][1]), 1 << BITS)
    weight, kappa = rational(coarse, 'original_far_weight'), rational(coarse, 'kappa')
    old_forward, old_lift = [rational(signed, key) for key in
                             ('sharp_forward_error_upper', 'lift_norm_upper')]
    new_lift = rational(coarse, 'new_direction_lift_norm_upper')
    hstar = rational(signed, 'moment_gram_norm_upper')
    if (kappa != fmpq(33, 32) or kappa != rational(signed, 'kappa')
            or weight != rational(signed, 'original_far_weight')
            or weight != fmpq(int(far['D_numerator']), 1 << far['D_grid_bits'])
            or hstar != fmpq(int(far['upper_bounds']['moment_gram_norm']), 1 << 628)
            or coarse['packet_sha256']['full-gamma-signed-correction.json.zlib'] !=
                PINS['full-gamma-signed-correction.json.zlib']):
        raise ValueError('Different original weight, uniform forward coefficient or moment trace.')
    width_ok = lambda pair: len(pair) == 2 and 0 <= int(pair[1])-int(pair[0]) <= 2
    if (len(far['moment_pairs']) != R or len(far['zeta']) != R
            or {row['U_power'] for row in far['moment_pairs']} != set(range(2, 2*R+1, 2))
            or any(row['T_power'] != row['U_power']+1 for row in far['moment_pairs'])
            or {row['p'] for row in far['zeta']} != set(range(4, 2*R+3, 2))
            or any(not width_ok(row[key]) for row in far['moment_pairs'] for key in ('U', 'T'))
            or any(not width_ok(row['interval']) for row in far['zeta'])
            or not width_ok(far['alpha_pair'])
            or len(coarse['new_short_jet_coefficients']) != 2*R
            or any(len(row) != 1 for row in coarse['new_short_jet_coefficients'])):
        raise ValueError('Missing saved new features or original midpoint/anchor contract.')
    # Metadata-only review stops before any actual feature/anchor conversion.
    forward = scaled_forward(old_forward, old_lift, new_lift)
    h = [interval_ball(row[0], BITS) for row in coarse['new_short_jet_coefficients']]
    anchor_row = next(row for row in far['moment_pairs'] if row['U_power'] == 2*R)
    anchor_zeta = next(row for row in far['zeta'] if row['p'] == 2*R+2)
    anchor = [interval_ball(pair, BITS) for pair in
              (anchor_row['U'], anchor_row['T'], anchor_zeta['interval'])]
    exact = directional_norm(h, upper(coarse, 'new_direction_gram_quadratic'),
        freeze_upper(interval_ball(primal['bounds']['actual_symbol_budget'], BITS).upper()),
        fmpq(int(far['alpha_pair'][1]), 1 << BITS),
        fmpq(int(far['alpha_pair'][1])-int(far['alpha_pair'][0]), 1 << (BITS+1)),
        weight, fmpq(2*R, 1 << BITS), hstar, forward, R, anchor, fmpq(M, M+1))
    old_norm2 = rational(coarse, 'new_direction_physical_norm2_upper')
    norm2 = min(old_norm2, exact['rebuilt_representer_energy_upper'])
    theta = ball(signed, 'actual_signed_theta')
    _, projected = projected_energy(fmpq(0), rational(coarse, 'source_energy_upper'),
        fmpq(0), norm2, rational(coarse, 'old_representer_norm2_upper'),
        theta/arb(kappa), ball(coarse, 'old_new_direction_physical_cross'))
    lo = max(rational(coarse, 'zero_pairing_lower'),
             freeze_lower(arb(kappa)*arb(projected['principal_mixed_cross_lower'])))
    hi = min(rational(coarse, 'zero_pairing_upper'),
             freeze_upper(arb(kappa)*arb(projected['principal_mixed_cross_upper'])))
    if lo > hi:
        raise ArithmeticError('Refined source interval contradicts its saved coarse interval.')
    zeta = arb(arb((lo+hi)/2), arb((hi-lo)/2))
    relaxed_zeta = ball(coarse, 'constructed_relaxation_zero_pairing')
    disjoint = bool(zeta.upper() < relaxed_zeta.lower() or relaxed_zeta.upper() < zeta.lower())
    derivative = zeta+ball(coarse, 'coordinate_ratio')*theta
    exact.update(new_sharp_forward_error_upper=forward, selected_direction_norm2_upper=norm2,
        saved_direction_norm2_upper=old_norm2, zero_pairing_lower=lo, zero_pairing_upper=hi,
        zero_pairing_absolute_upper=max(abs(lo), abs(hi)), original_far_weight=weight, kappa=kappa,
        source_energy_upper=rational(coarse, 'source_energy_upper'))
    scalars = {key: arb(value) for key, value in exact.items()}
    scalars.update(actual_far_zero_pairing=zeta,
        actual_far_reference_rho_pairing=derivative,
        constructed_relaxation_zero_pairing=relaxed_zeta)
    require_finite(list(scalars.values()))
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    if any(pair is None or int(pair[0]) > int(pair[1]) for pair in exported.values()):
        raise ArithmeticError('Invalid directional refinement export; stop without retry.')
    result = dict(schema='ACTUAL_C27_ZERO_CHANNEL_JET_REFINEMENT_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: digest for p, digest in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: digest for p, digest in guards.items() if p not in invoked},
        source_integer_sha256=coarse['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=1024, scalar_bounds=exported,
        exact_rationals={key: [str(value.p), str(value.q)] for key, value in exact.items()},
        norm_strictly_improved=bool(norm2 < old_norm2),
        source_pairing_strictly_narrower=bool(max(abs(lo), abs(hi)) <
            rational(coarse, 'zero_pairing_absolute_upper')),
        constructed_relaxation_row_excluded_by_source_pairing=disjoint,
        new_directional_norm_applications=1, new_homogeneous_forward_bound_applications=1,
        new_projection_applications=1, new_direction_feature_quadratic_or_cross_actions=0,
        new_source_sine_moment_tail_inverse_solve_or_trace_suppliers=0, old_producers_replayed=0,
        parameter_source_precision_cutoff_scan=False,
        actual_center_sign_claim=False, actual_full_border_sign_claim=False, actual_rh_proof=False,
        conditions=coarse['conditions']+[
            'The original sharp weighted forward allowance equals a uniform direction-independent coefficient times the original positive lift cap Lb, up to outward rounding. Thus old_forward*Ld/Lb is valid for the SAME new direction; no arbitrary realization-specific norm cap is rescaled.',
            'The existing directional_norm applies to the SAVED new features and quadratic on the SAME moment mask. Three retained anchors, omission and rounding payments are reused as a method, not old numerical outputs or new moment suppliers.',
            'The saved old/new physical cross remains its existing interval. Its intersection with the new squared-norm cap is valid for the SAME source and directions. No old direction, feature, kernel quadratic or cross is regenerated.',
            'A disjoint interval excludes only the saved constructed relaxation row under the added source relation, not all negative rows, a physically legal alpha, the actual paid comparison, Robin or RH.'
        ], scope='One saved-output application of existing directional moment, homogeneous forward and common-projection bounds; no generic theorem, novelty or Lean claim.',
        elapsed_seconds=time.monotonic()-started)
    verify_guards(guards)
    digest = write_packet(out, result)
    print('saved', out.name, digest, flush=True)
    for key in ('selected_direction_norm2_upper', 'new_sharp_forward_error_upper',
                'actual_far_zero_pairing', 'actual_far_reference_rho_pairing',
                'constructed_relaxation_zero_pairing'):
        print(key, scalars[key].str(24), flush=True)
    print('narrower', result['source_pairing_strictly_narrower'],
          'constructed row excluded', disjoint, 'seconds', time.monotonic()-started, flush=True)


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
