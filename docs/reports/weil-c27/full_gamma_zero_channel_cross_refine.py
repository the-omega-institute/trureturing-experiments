"""Recover a physical cross interval from its saved masked contraction.

Use two-sided directional moment and forward payments on saved features.
Apply existing common projection and scalar square completion; no finite
direction, feature, quadratic, cross or original supplier is repeated.
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
from full_gamma_negative_index_joint import freeze_lower
from full_gamma_shared_projection_energy import projected_energy
from full_gamma_signed_mu import require_finite, verify_guards, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS, R, M = 640, 56, 32768
PINS = {
    'full-gamma-zero-channel-jet-refine.json.zlib':
        'ee401d13ef503a44089c33e958b8f90dab09aefa843050f09713ee6a33729861',
    'full-gamma-zero-channel-pairing.json.zlib':
        'c4cbbead777f652812c0aa0bd238bf08053daa52be7afbedfa1321da19847a9b',
    'full-gamma-directional-jet-norm.json.zlib':
        '9ccd6494c0520f9f53efe0b44bcae879566157c7e7ffb698165f9d02c258fd41',
    'full-gamma-signed-jet-interface.json.zlib':
        '154949d8ab3e1f2875de531aa63a9c5a3c3d7d133640a5ce8b1d23c98d67565a',
    'full-gamma-signed-correction.json.zlib':
        'd527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd',
    'far-coupling.json.zlib':
        'ac35edacd7f55de470586e3b8a5915bcaa49b9d2c4c9c2ee073eeaaa8c9c73e7',
    'full-gamma-primal-tail.json':
        'fad634fb25d8d20e83a15d8277569978646246c2a560d0dc0f7a854b5087455a',
    'full-gamma-steering-capacity.json.zlib':
        '85a52adf97e96def95aaf5abaa307338522a2bae5d7b40c96692ad0366a444fe',
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
}


def bilinear_payment(hb, hd, symbol, weight, alpha, radius, moment_error,
                     trace, hb2, hd2, jb2, jd2, eb, ed, diagonals):
    """Pay BOTH signs of the omitted actual-symbol jet cross."""
    order = len(hb)//2
    vals = (symbol, weight, alpha, radius, moment_error, trace,
            hb2, hd2, jb2, jd2, eb, ed)
    require_finite(list(hb)+list(hd)+[arb(x) for x in vals])
    if (order < 1 or len(hb) != 2*order or len(hd) != len(hb)
            or type(diagonals) is not int or not 1 <= diagonals <= 2*order-1
            or not symbol >= 1 or not weight > 0 or not alpha > 0
            or any(not v >= 0 for v in vals[3:])):
        raise ValueError('Missing same-jet shapes or nonnegative bilinear caps.')
    total = fmpq(0)
    for n in range(diagonals, 2*order-1):
        pairs = [(i, n-i) for i in range(order) if 0 <= n-i < order]
        aa = sum((hb[i]*hd[j] for i, j in pairs), arb(0))
        ab = sum((hb[i]*hd[order+j]+hb[order+i]*hd[j]
                  for i, j in pairs), arb(0))
        bb = sum((hb[order+i]*hd[order+j] for i, j in pairs), arb(0))
        total += freeze_upper(abs(aa).upper()*arb(symbol)**2/(2*n+1)
            +abs(ab).upper()*arb(symbol)/(2*n+2)+abs(bb).upper()/(2*n+3))
    masking = freeze_upper(arb(weight)*arb(alpha)*arb(total))
    rounding = freeze_upper(arb(weight)*(arb(hb2)*arb(hd2)).sqrt()
        *(arb(alpha)*arb(moment_error)+arb(radius)*(arb(trace)+arb(moment_error))))
    forward = freeze_upper(arb(eb)*arb(jd2).sqrt()+arb(ed)*arb(jb2).sqrt()
                           +arb(eb)*arb(ed))
    return dict(bilinear_omitted_moment_absolute_upper=total,
        bilinear_masking_error_upper=masking, bilinear_rounding_error_upper=rounding,
        bilinear_forward_error_upper=forward,
        complete_bilinear_cross_error_upper=masking+rounding+forward)


def paid_floor(alpha0, k0, fee0, q0, d0, h11, theta, response,
               negative_credit, constant_loss, curvature_loss, k, fee):
    """Reuse saved comparison losses with a new signed response interval."""
    vals = (alpha0, k0, fee0, q0, d0, h11, theta, response,
            negative_credit, constant_loss, curvature_loss, k, fee)
    require_finite([arb(x) for x in vals])
    if (min(k0, fee0, k, fee) <= 0
            or min(negative_credit, constant_loss, curvature_loss) < 0):
        raise ValueError('Missing old/new pricing or saved loss domains.')
    shift = k-k0
    base = freeze_lower(q0.lower()+arb(fee-fee0)+arb(shift)*arb(alpha0)**2
        +2*theta.lower()-arb(negative_credit)-arb(constant_loss))
    curvature = freeze_lower(h11.lower()+arb(shift)-arb(curvature_loss))
    slope = freeze_upper(abs(d0+arb(shift)*arb(alpha0)+response).upper())
    result = dict(comparison_base_lower=base, comparison_curvature_lower=curvature,
                  total_slope_absolute_upper=slope)
    if curvature > 0:
        gain = freeze_upper(arb(slope)**2/arb(curvature))
        result.update(maximal_steering_gain_upper=gain, all_alpha_floor_lower=base-gain)
    return result


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
            raise ValueError('Different accepted bilinear input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    new, coarse, latest, jet, signed, far, primal, steering, old = packets.values()
    inherited = dict(PINS)
    for packet in (new, coarse, latest, steering):
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256'):
            for name, digest in packet[field].items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting bilinear dependency: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked bilinear dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_ZERO_CHANNEL_JET_REFINEMENT_V1',
            'ACTUAL_C27_ZERO_CHANNEL_PAIRING_BOUND_V1', 'ACTUAL_C27_DIRECTIONAL_JET_NORM_V1',
            'ACTUAL_C27_SIGNED_JET_INTERFACE_V1', 'ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1',
            'ACTUAL_C27_FAR_COUPLING_GRAM_V1', 'ACTUAL_C27_FULL_GAMMA_PRIMAL_TAIL_V1',
            'ACTUAL_C27_ALL_ALPHA_COMPARISON_CAPACITY_V1',
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any(p['source_integer_sha256'] != new['source_integer_sha256']
                   or (p['source_bits'], p['source_row']) != (384, 0)
                   for p in (new, coarse, latest, jet, signed, steering, old))
            or (far['c'], far['N'], far['M'], far['jet_order'],
                far['moment_diagonals'], far['bound_grid_bits']) != (27, 1280, M, R, R, 628)
            or (jet['physical_last_mode'], jet['first_far_mode'], jet['jet_order'],
                jet['feature_order'], jet['original_metadata_grid_bits']) !=
                (2048, M+1, R, 'A0..A55,B0..B55', 512)
            or primal['metadata']['input_sha256'] != jet['input_sha256']
            or primal['metadata']['basis'] !=
                '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'):
        raise ValueError('Different same-source physical/moment/reference interface.')
    rational = lambda p, key: fmpq(*map(int, p['exact_rationals'][key]))
    ball = lambda p, key: interval_ball(p['scalar_bounds'][key], BITS)
    lower = lambda p, key: fmpq(int(p['scalar_bounds'][key][0]), 1 << BITS)
    upper = lambda p, key: fmpq(int(p['scalar_bounds'][key][1]), 1 << BITS)
    weight, kappa = rational(new, 'original_far_weight'), rational(new, 'kappa')
    alpha, radius = [rational(latest, key) for key in ('alpha_upper', 'alpha_midpoint_radius')]
    moment_error = rational(latest, 'moment_midpoint_error_upper')
    hstar = rational(signed, 'moment_gram_norm_upper')
    if (kappa != fmpq(33, 32) or kappa != rational(signed, 'kappa')
            or weight != rational(signed, 'original_far_weight')
            or weight != fmpq(int(far['D_numerator']), 1 << far['D_grid_bits'])
            or alpha != fmpq(int(far['alpha_pair'][1]), 1 << BITS)
            or radius != fmpq(int(far['alpha_pair'][1])-int(far['alpha_pair'][0]), 1 << (BITS+1))
            or moment_error != fmpq(2*R, 1 << BITS)
            or hstar != fmpq(int(far['upper_bounds']['moment_gram_norm']), 1 << 628)
            or rational(new, 'source_energy_upper') != rational(latest, 'selected_far_energy_upper')
            or rational(coarse, 'old_representer_norm2_upper') !=
                rational(latest, 'selected_representer_energy_upper')
            or any(len(p[key]) != 2*R or any(len(row) != 1 for row in p[key])
                   for p, key in ((jet, 'short_jet_coefficients'),
                                  (coarse, 'new_short_jet_coefficients')))
            or rational(steering, 'alpha0') != rational(old, 'fixed_alpha')
            or rational(steering, 'legal_coefficient') != rational(old, 'legal_coefficient')
            or rational(steering, 'frozen_positive_fee') != rational(old, 'positive_source_valid_fee')):
        raise ValueError('Different saved features, moments, source energy or pricing origin.')
    # Metadata-only source review stops before any actual feature conversion.
    hb = [interval_ball(row[0], BITS) for row in jet['short_jet_coefficients']]
    hd = [interval_ball(row[0], BITS) for row in coarse['new_short_jet_coefficients']]
    exact = bilinear_payment(hb, hd,
        freeze_upper(interval_ball(primal['bounds']['actual_symbol_budget'], BITS).upper()),
        weight, alpha, radius, moment_error, hstar,
        rational(latest, 'feature_absolute_norm_squared_upper'),
        rational(new, 'feature_absolute_norm_squared_upper'),
        rational(latest, 'full_jet_energy_upper'), rational(new, 'full_jet_energy_upper'),
        rational(signed, 'sharp_forward_error_upper'),
        rational(new, 'new_sharp_forward_error_upper'), R)
    error = exact['complete_bilinear_cross_error_upper']
    masked = ball(coarse, 'old_new_direction_kernel_cross')
    clo = max(lower(coarse, 'old_new_direction_physical_cross'), freeze_lower(masked.lower()-arb(error)))
    chi = min(upper(coarse, 'old_new_direction_physical_cross'), freeze_upper(masked.upper()+arb(error)))
    if clo > chi:
        raise ArithmeticError('Bilinear certificate contradicts the saved cross interval.')
    cross = arb(arb((clo+chi)/2), arb((chi-clo)/2))
    energy = rational(new, 'source_energy_upper')
    jb, jd = rational(coarse, 'old_representer_norm2_upper'), rational(new, 'selected_direction_norm2_upper')
    theta = ball(signed, 'actual_signed_theta')
    _, projection = projected_energy(fmpq(0), energy, fmpq(0), jd, jb, theta/arb(kappa), cross)
    zlo = max(rational(new, 'zero_pairing_lower'),
              freeze_lower(arb(kappa)*arb(projection['principal_mixed_cross_lower'])))
    zhi = min(rational(new, 'zero_pairing_upper'),
              freeze_upper(arb(kappa)*arb(projection['principal_mixed_cross_upper'])))
    if zlo > zhi:
        raise ArithmeticError('Bilinear projection contradicts the saved source interval.')
    zeta = arb(arb((zlo+zhi)/2), arb((zhi-zlo)/2))
    response = zeta+ball(coarse, 'coordinate_ratio')*theta
    rlo = max(lower(new, 'actual_far_reference_rho_pairing'), freeze_lower(response.lower()))
    rhi = min(upper(new, 'actual_far_reference_rho_pairing'), freeze_upper(response.upper()))
    if rlo > rhi:
        raise ArithmeticError('Bilinear response contradicts the saved reference pairing.')
    response = arb(arb((rlo+rhi)/2), arb((rhi-rlo)/2))
    alpha0, k0, fee0 = (rational(steering, 'alpha0'), 1/rational(steering, 'legal_coefficient'),
                        rational(steering, 'frozen_positive_fee'))
    for prefix, k, fee in (
            ('fixed', 1/rational(latest, 'repriced_legal_coefficient'),
             rational(latest, 'source_valid_fee_upper')),
            ('energy_endpoint', rational(latest, 'endpoint_reciprocal_coefficient_lower'),
             rational(latest, 'endpoint_comparison_fee'))):
        bounds = paid_floor(alpha0, k0, fee0, ball(steering, 'reference_finite_expression'),
            ball(steering, 'reference_finite_derivative'), ball(old, 'steering_H11_reference'),
            theta, response, rational(steering, 'negative_far_square_credit_upper'),
            rational(steering, 'uniform_comparison_constant_loss_upper'),
            rational(steering, 'uniform_comparison_curvature_loss_upper'), k, fee)
        exact.update({prefix+'_'+key: value for key, value in bounds.items()})
        exact.update({prefix+'_reciprocal_coefficient': k, prefix+'_comparison_fee': fee})
    exact.update(physical_cross_lower=clo, physical_cross_upper=chi,
        zero_pairing_lower=zlo, zero_pairing_upper=zhi,
        reference_rho_pairing_lower=rlo, reference_rho_pairing_upper=rhi,
        source_energy_upper=energy, old_representer_norm2_upper=jb,
        new_representer_norm2_upper=jd, original_far_weight=weight, kappa=kappa)
    scalars = {key: arb(value) for key, value in exact.items()}
    scalars.update(actual_physical_cross=cross, actual_far_zero_pairing=zeta,
                   actual_far_reference_rho_pairing=response)
    require_finite(list(scalars.values()))
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    if any(pair is None or int(pair[0]) > int(pair[1]) for pair in exported.values()):
        raise ArithmeticError('Invalid bilinear result export; stop without retry.')
    result = dict(schema='ACTUAL_C27_ZERO_CHANNEL_CROSS_REFINEMENT_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: digest for p, digest in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: digest for p, digest in guards.items() if p not in invoked},
        source_integer_sha256=new['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=1024, scalar_bounds=exported,
        exact_rationals={key: [str(value.p), str(value.q)] for key, value in exact.items()},
        physical_cross_strictly_narrower=bool(clo > lower(coarse, 'old_new_direction_physical_cross')
                                             or chi < upper(coarse, 'old_new_direction_physical_cross')),
        source_pairing_strictly_narrower=bool(zlo > rational(new, 'zero_pairing_lower')
                                             or zhi < rational(new, 'zero_pairing_upper')),
        fixed_paid_recipe_excluded=bool(exact.get('fixed_all_alpha_floor_lower', fmpq(0)) > 0),
        energy_cap_only_family_excluded=bool(exact.get('energy_endpoint_all_alpha_floor_lower', fmpq(0)) > 0),
        new_bilinear_error_applications=1, new_projection_applications=1,
        new_scalar_paid_floor_applications=2, new_finite_direction_feature_quadratic_or_cross_actions=0,
        new_source_sine_moment_tail_inverse_solve_or_trace_suppliers=0,
        old_producers_replayed=0, parameter_source_precision_cutoff_scan=False,
        actual_center_sign_claim=False, actual_full_border_sign_claim=False, actual_rh_proof=False,
        conditions=[
            'SAME actual-symbol full jets, physical source, reference, weights, original mask and forward certificates. No approximate-symbol finite samples are used.',
            'The saved masked cross is retained; absolute ordered AA/AB+BA/BB antidiagonals pay BOTH omitted signs, plus separate moment/prefactor and forward-vector errors. Rounded or masked kernels are not assumed PSD.',
            'Existing common projection retains the entire actual theta interval and SAME source energy E with directional norm uppers Jb,Jd. Empty intersections stop; no independently favorable extrema are asserted attainable.',
            'Existing square completion uses the saved old reference/negative credit/inverse-loss uppers, with explicit fee and reciprocal-coefficient translation. It covers all real alpha without asserting physical legality of every alpha.',
            'The saved energy LOWER-endpoint kmin,Dmin bound the fixed monotone energy-cap-only pricing rule. A positive floor excludes only that recipe family; a nonpositive lower is inconclusive, not a realized negative row or actual-source upper.',
            'Actual-center/full-border/good-block/cofinal positivity, Robin and RH remain unresolved. No generic theorem, novelty or Lean proof is claimed.'
        ], elapsed_seconds=time.monotonic()-started)
    verify_guards(guards)
    digest = write_packet(out, result)
    print('saved', out.name, digest, flush=True)
    for key in ('complete_bilinear_cross_error_upper', 'actual_physical_cross',
                'actual_far_zero_pairing', 'actual_far_reference_rho_pairing',
                'fixed_all_alpha_floor_lower', 'energy_endpoint_all_alpha_floor_lower'):
        if key in scalars:
            print(key, scalars[key].str(24), flush=True)
    print('narrower', result['physical_cross_strictly_narrower'],
          result['source_pairing_strictly_narrower'], 'recipes excluded',
          result['fixed_paid_recipe_excluded'], result['energy_cap_only_family_excluded'],
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
