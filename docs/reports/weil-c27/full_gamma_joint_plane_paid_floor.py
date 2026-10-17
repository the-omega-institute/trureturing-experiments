"""Retain two saved signed responses in a negative-index-one projection.

Apply the classical two-dimensional Schur bound and a clipped scalar
quadratic minimum. Saved scalars only; no original supplier is repeated.
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
from full_gamma_signed_mu import require_finite, verify_guards, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
PINS = {
    'full-gamma-zero-channel-cross-refine.json.zlib':
        '5cfba16d5428c1c61618a861a8167ca6828f46142d7bd7993d4d84b20de6e980',
    'full-gamma-zero-channel-pairing.json.zlib':
        'c4cbbead777f652812c0aa0bd238bf08053daa52be7afbedfa1321da19847a9b',
    'full-gamma-relaxation-row-witness.json.zlib':
        '1879f4b0653e762ca0d7bb1e1933a59afc79fb10851d155f1c16c87055872992',
    'full-gamma-negative-index-joint.json.zlib':
        'a7df6338b629e9dd791ac2914c8efc7f3ef3c69c20e562da0aa18067aa40728d',
    'full-gamma-signed-correction.json.zlib':
        'd527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd',
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
    'full-gamma-directional-jet-norm.json.zlib':
        '9ccd6494c0520f9f53efe0b44bcae879566157c7e7ffb698165f9d02c258fd41',
}


def clipped_quadratic(a, linear, constant, zlo, zhi):
    """Lower an interval quadratic on its entire closed rational interval."""
    require_finite([arb(v) for v in (a, linear, constant, zlo, zhi)])
    if not a > 0 or zlo > zhi:
        raise ValueError('A positive rational curvature and ordered interval are required.')
    blo, bhi = freeze_lower(linear.lower()), freeze_upper(linear.upper())
    midpoint, radius = (blo+bhi)/2, (bhi-blo)/2
    root = max(zlo, min(zhi, -midpoint/a))
    debit = 2*radius*max(abs(zlo), abs(zhi))
    base = freeze_lower(constant.lower())
    floor = base+a*root*root+2*midpoint*root-debit
    return dict(quadratic_coefficient_lower=a, linear_midpoint=midpoint,
                linear_radius_upper=radius, constant_lower=base,
                interval_linear_error_upper=debit, clipped_model_zeta=root,
                box_floor_lower=floor)


def plane_floor(alpha0, b, g, h, c, t, theta, zlo, zhi, gamma, k, fee):
    """Project the SAME form onto y0,rho, then lower over both response boxes.

    The inputs have one common realization with b=-y0 Qeps y0*>0,
    g=rho Qeps y0*, h=rho Qeps rho*, t=w Qeps w*>=0 and
    w=rho-c*y0. These identities are inherited premises, not inferred
    from overlapping intervals. Only scalar interval arithmetic is new.
    """
    require_finite([arb(v) for v in
        (alpha0, b, g, h, c, t, theta, zlo, zhi, gamma, k, fee)])
    if (not b > 0 or not t >= 0 or not 0 < gamma <= 1
            or not k > 0 or not fee > 0 or zlo > zhi):
        raise ValueError('Missing common timelike/positive/pricing/interval domains.')
    omega = c+g/b
    a = t+b*omega*omega
    same_a = h+g*g/b
    if not a.overlaps(same_a):
        raise ArithmeticError('Inconsistent saved common-plane identities.')
    curvature, rank_gap = arb(k)+h, arb(k)-g*g/b
    values = dict(positive_plane_energy=a, comparison_curvature=curvature,
                  rank_gap=rank_gap, transported_theta=theta)
    exact = {}
    if not a > 0 or not curvature > 0 or not rank_gap > 0:
        return values, exact
    coefficient = arb(gamma)**2*(1/a-1/curvature)
    a_lower = freeze_lower(coefficient.lower())
    if not a_lower > 0:
        return values, exact
    theta_lo, theta_hi = freeze_lower(theta.lower()), freeze_upper(theta.upper())
    exact.update(theta_lower=theta_lo, theta_upper=theta_hi)
    floors = []
    for prefix, endpoint in (('theta_lo', theta_lo), ('theta_hi', theta_hi)):
        th = arb(endpoint)
        slope = arb(k)*arb(alpha0)+g+c*th
        linear = arb(gamma)*(omega*th/a-slope/curvature)
        constant = arb(k)*arb(alpha0)**2+arb(fee)-b+2*th-th*th/b \
            +(omega*th)**2/a-slope*slope/curvature
        bound = clipped_quadratic(a_lower, linear, constant, zlo, zhi)
        exact.update({prefix+'_'+key: value for key, value in bound.items()})
        floors.append(bound['box_floor_lower'])
    exact['all_alpha_floor_lower'] = min(floors)
    return values, exact


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
            raise ValueError('Different accepted joint-plane input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw))
    cross, coarse, relaxed, joint, signed, old, latest = [packets[n] for n in PINS]
    inherited = dict(PINS)
    for packet in (cross, coarse, relaxed, joint, latest):
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256'):
            for name, digest in packet[field].items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting joint-plane dependency: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked joint-plane dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_ZERO_CHANNEL_CROSS_REFINEMENT_V1',
            'ACTUAL_C27_ZERO_CHANNEL_PAIRING_BOUND_V1',
            'C27_SAVED_CONSTRAINT_RELAXATION_ROW_V1',
            'ACTUAL_C27_NEGATIVE_INDEX_JOINT_BOUND_V1',
            'ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1',
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1',
            'ACTUAL_C27_DIRECTIONAL_JET_NORM_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any(p['source_integer_sha256'] != signed['source_integer_sha256']
                   or (p['source_bits'], p['source_row']) != (384, 0)
                   for p in (cross, coarse, joint, old, latest))
            or relaxed['reference_source_integer_sha256'] != signed['source_integer_sha256']
            or (relaxed['reference_source_bits'], relaxed['reference_source_row']) != (384, 0)):
        raise ValueError('Different SAME-source scalar roles or grids.')
    rational = lambda p, key: fmpq(*map(int, p['exact_rationals'][key]))
    upper = lambda p, key: fmpq(int(p['scalar_bounds'][key][1]), 1 << BITS)
    if (rational(relaxed, 'reference_alpha0') != rational(old, 'fixed_alpha')
            or rational(relaxed, 'source_energy_upper') != rational(cross, 'source_energy_upper')
            or rational(cross, 'fixed_reciprocal_coefficient') !=
                1/rational(latest, 'repriced_legal_coefficient')
            or rational(cross, 'fixed_comparison_fee') != rational(latest, 'source_valid_fee_upper')
            or rational(cross, 'energy_endpoint_reciprocal_coefficient') !=
                rational(latest, 'endpoint_reciprocal_coefficient_lower')
            or rational(cross, 'energy_endpoint_comparison_fee') !=
                rational(latest, 'endpoint_comparison_fee')
            or rational(cross, 'kappa') != fmpq(33, 32)
            or rational(cross, 'kappa') != rational(signed, 'kappa')):
        raise ValueError('Different reference, source energy or saved paid recipes.')
    # Metadata-only source review stops before this actual scalar conversion.
    ball = lambda p, key: interval_ball(p['scalar_bounds'][key], BITS)
    alpha0 = rational(relaxed, 'reference_alpha0')
    eta = upper(old, 'neumann_eta_upper')
    if not 0 <= eta < 1:
        raise ValueError('Missing accepted two-sided inverse domain.')
    epsilon, gamma = eta/(1+eta), 1/(1+eta)
    r0 = ball(old, 'selected_reference_form')
    rr = ball(old, 'steering_reference_rRr')
    p = ball(relaxed, 'reference_rho_positive_pairing')
    b = ball(relaxed, 'reference_timelike_energy')+arb(epsilon)*r0
    g = ball(joint, 'reference_rho_pairing')-arb(epsilon)*p
    h = ball(old, 'steering_H11_reference')-arb(1/rational(old, 'legal_coefficient')) \
        -arb(epsilon)*rr
    c = ball(coarse, 'coordinate_ratio')
    t = arb(gamma)*ball(coarse, 'zero_direction_positive_energy')
    afar = rational(relaxed, 'relaxed_positive_energy_cap')
    if not r0 >= 0 or not afar >= 0:
        raise ValueError('Missing positive-form caps for the SAME actual far row.')
    theta_radius = freeze_upper(arb(epsilon)*(arb(afar)*r0.upper()).sqrt().upper())
    theta = ball(signed, 'actual_signed_theta')+arb(0, arb(theta_radius))
    zlo, zhi = rational(cross, 'zero_pairing_lower'), rational(cross, 'zero_pairing_upper')
    scalars = dict(timelike_energy=b, rho_timelike_pairing=g,
        rho_reference_energy=h, coordinate_ratio=c, zero_coordinate_energy=t)
    exact = dict(epsilon=epsilon, gamma=gamma, alpha0=alpha0,
        transported_theta_error_upper=theta_radius, fixed_far_positive_energy_upper=afar,
        zeta_lower=zlo, zeta_upper=zhi)
    for prefix in ('fixed', 'energy_endpoint'):
        k, fee = (rational(cross, prefix+'_'+key) for key in
                  ('reciprocal_coefficient', 'comparison_fee'))
        values, bounds = plane_floor(alpha0, b, g, h, c, t, theta, zlo, zhi, gamma, k, fee)
        scalars.update({prefix+'_'+key: value for key, value in values.items()})
        exact.update({prefix+'_'+key: value for key, value in bounds.items()})
        previous = rational(latest, prefix+'_all_alpha_floor_lower')
        retained = max(previous, bounds.get('all_alpha_floor_lower', previous))
        exact.update({prefix+'_retained_all_alpha_floor_lower': retained,
                      prefix+'_previous_all_alpha_floor_lower': previous,
                      prefix+'_reciprocal_coefficient': k, prefix+'_comparison_fee': fee})
    scalars.update({key: arb(value) for key, value in exact.items()})
    require_finite(list(scalars.values()))
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    if any(pair is None or int(pair[0]) > int(pair[1]) for pair in exported.values()):
        raise ArithmeticError('Invalid joint-plane export; stop without retry.')
    improved = lambda p: exact[p+'_retained_all_alpha_floor_lower'] > exact[p+'_previous_all_alpha_floor_lower']
    excluded = lambda p: exact[p+'_retained_all_alpha_floor_lower'] > 0
    result = dict(schema='ACTUAL_C27_JOINT_PLANE_PAID_FLOOR_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: digest for p, digest in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: digest for p, digest in guards.items() if p not in invoked},
        source_integer_sha256=signed['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=1024, scalar_bounds=exported,
        exact_rationals={key: [str(value.p), str(value.q)] for key, value in exact.items()},
        fixed_floor_strictly_improved=bool(improved('fixed')),
        endpoint_floor_strictly_improved=bool(improved('energy_endpoint')),
        fixed_paid_recipe_excluded=bool(excluded('fixed')),
        energy_cap_only_family_excluded=bool(excluded('energy_endpoint')),
        new_joint_plane_applications=2, theta_endpoints_per_recipe=2,
        old_producers_or_fixtures_replayed=0, new_finite_or_physical_suppliers=0,
        parameter_source_precision_cutoff_scan=False,
        actual_center_sign_claim=False, actual_full_border_sign_claim=False, actual_rh_proof=False,
        conditions=[
            'All finite scalars and both signed responses have the SAME original source/core/lift/reference/real-form interpretation. Original actual-symbol jet/moment/forward certificates remain premises.',
            'Qeps=(1-epsilon)P-e0e0*/delta <= Qtrue, P>=0; y0 is timelike and w=rho-c*y0 has zero first coordinate. The two-dimensional Schur projection retains its positive square, including singular P.',
            'The theta transport uses a saved UPPER on the SAME actual fPf*, never a source-energy lower. Source-response boxes are outer relaxations, not claims of joint endpoint attainability.',
            'Concavity in transported theta reduces its minimum to two endpoints. A positive rank gap makes the zeta quadratic convex. Directed coefficient bounds and a linear-radius debit give an exact clipped model minimum; its optimizer need not be an actual row or physical steering value.',
            'Both recipes retain their saved fee and coefficient. A positive endpoint floor excludes only the fixed monotone energy-cap-only pricing family, not a smaller-charge actual center. A nonpositive lower supplies no negative realized value or upper.',
            'Stronger previous floors remain retained. Actual-center/full-border/good-block/cofinal positivity, Robin and RH remain unresolved. No generic theorem, novelty or Lean proof is claimed.'
        ], elapsed_seconds=time.monotonic()-started)
    verify_guards(guards)
    digest = write_packet(out, result)
    print('saved', out.name, digest, flush=True)
    for prefix in ('fixed', 'energy_endpoint'):
        key = prefix+'_retained_all_alpha_floor_lower'
        print(key, scalars[key].str(24), 'improved', bool(improved(prefix)),
              'excluded', bool(excluded(prefix)), flush=True)


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
