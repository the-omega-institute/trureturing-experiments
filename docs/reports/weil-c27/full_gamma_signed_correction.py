"""Supply the missing SAME-source mixed correction from the saved Gram kernel.

Reuse the completed signed proxy, full112 jet, physical Gram certificate,
inverse profiles and separate forward bounds. No old producer is replayed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, arb_mat, ctx, fmpq

from arb_inverse_residual_obstruction import inputs
from diagonal_floor_witness import interval_ball
from full_gamma_coherent_split_capacity import freeze_upper
from full_gamma_signed_mu import (require_finite, source_polynomial,
                                  verify_guards, write_packet)
from joint_gamma_boundary import matrix_in, matrix_out
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS, M, N, R, FIRST = 640, 32768, 2048, 56, 1281
PINS = {
    'full-gamma-signed-mu-direct-sine.json.zlib':
        '6475e5f63dba423b2f1ee7e488581f15ef8d53de83a22ea064ac2e393b124d55',
    'full-gamma-signed-jet-interface.json.zlib':
        '154949d8ab3e1f2875de531aa63a9c5a3c3d7d133640a5ce8b1d23c98d67565a',
    'structured-joint-2048-target.json.zlib':
        'ea4a761f76b59a11b2d04faf80531901a5164f4145c2424a4e79889ae3f72e22',
    'far-coupling.json.zlib':
        'ac35edacd7f55de470586e3b8a5915bcaa49b9d2c4c9c2ee073eeaaa8c9c73e7',
    'full-gamma-dual-32768-defect.json.zlib':
        '6f320576338562a4e7338906ca66b25ab49e4a1d4460e87e26165cd9579a1154',
    'full-gamma-primal-tail.json':
        'fad634fb25d8d20e83a15d8277569978646246c2a560d0dc0f7a854b5087455a',
    'full-gamma-hankel-tail-energy.json.zlib':
        'd09a5ab254d164cfd40982b31e1a155e31d6288b705d916bf248236a0d527606',
}


def left_features(values, midpoints, first, scale, order):
    """Apply the existing signed power reducer to the two P0 feature families."""
    if len(values) != len(midpoints) or not values:
        raise ValueError('Require the same nonempty source and midpoint rows.')
    signed = [(-1)**(first+i)*v for i, v in enumerate(values)]
    apoly = source_polynomial(signed, first, scale, 2*order)
    bpoly = source_polynomial([-v*s for v, s in zip(signed, midpoints)],
                              first, scale, 2*order)
    return arb_mat([[apoly[2*k]] for k in range(order)]
                   +[[bpoly[2*k+1]] for k in range(order)])


def price_correction(b_lower, b_upper, weight, hmix, kphys, lift_norm,
                     eps_far, hstar, e0, rho, eps_sm, ncap, principal, kappa,
                     physical_last=N, scale=M, order=R):
    """Apply the inherited mixed norm formula and separate forward components."""
    values = [arb(v) for v in (b_lower, b_upper, weight, hmix, kphys, lift_norm,
                              eps_far, hstar, e0, rho, eps_sm, ncap, principal, kappa)]
    require_finite(values)
    if (not b_lower > 0 or not b_upper >= b_lower or not weight > 0
            or not kphys > 0 or not kappa > 0
            or any(not v >= 0 for v in values[3:-1])
            or not 0 < physical_last < scale or order < 1):
        raise ValueError('Invalid positive source/coupling/forward bound domain.')
    pi, ratio = arb.pi(), arb(fmpq(physical_last, scale))
    c_lower = (arb(weight).sqrt()*2*b_lower
               /(pi*arb(scale-physical_last).sqrt())).lower()
    if not c_lower > 0:
        raise ArithmeticError('The original mixed-norm coefficient is not positive.')
    # Hmix >= c_actual*||lambda||_1 and c_actual >= c_lower > 0.
    z_lambda = freeze_upper(arb(hmix)/c_lower)
    separation = 1-arb(fmpq(physical_last, scale+1))
    r2 = 4*b_upper**2*physical_last/(pi**2*separation**2*scale*(4*order+1))
    r2 *= ratio**(4*order)
    ep2 = arb(e0)**2+arb(rho)**2*physical_last*ratio**2/(1-ratio**4)
    alpha = 4/(pi**2*scale)
    # r2 is below the640 grid: take its sqrt at working precision first.
    eps_b = freeze_upper((arb(weight)*arb(kphys)).sqrt()
        *(r2.sqrt()+(alpha*arb(hstar)*ep2).sqrt())*arb(lift_norm))
    cross_error = freeze_upper(arb(kappa)*arb(eps_far)*arb(z_lambda)
                              *arb(kphys).sqrt()*arb(lift_norm))
    smooth_error = freeze_upper(arb(kappa)*arb(eps_sm)*arb(ncap))
    forward_error = freeze_upper(arb(kappa)*arb(principal)*arb(eps_b))
    total = cross_error+smooth_error+forward_error
    exact = dict(lambda_l1_upper=z_lambda, sharp_forward_error_upper=eps_b,
                 mixed_cross_error_upper=cross_error, smooth_error_upper=smooth_error,
                 proxy_forward_error_upper=forward_error, complete_theta_error_upper=total)
    scalars = {k: arb(v) for k, v in exact.items()}
    scalars.update(mixed_norm_coefficient_lower=c_lower,
                   feature_frobenius_error_upper=ep2.sqrt(),
                   forward_jet_remainder_norm_upper=r2.sqrt())
    require_finite(list(scalars.values())+[r2])
    return exact, scalars, r2


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
            raise ValueError('Different accepted correction packet: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    mu, jet, target, far, defect, primal, hankel = [packets[n] for n in PINS]
    inherited = dict(PINS)
    for packet in (mu, jet):
        for field in ('packet_sha256', 'invocation_source_sha256', 'inherited_dependency_sha256'):
            for name, digest in packet[field].items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting correction input: '+name)
                inherited[name] = digest
    for name, digest in target['metadata']['input_sha256'].items():
        if name in inherited and inherited[name] != digest:
            raise ValueError('Conflicting physical scalar input: '+name)
        inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if [p['schema'] for p in (mu, jet, target, far, defect, primal, hankel)] != [
            'ACTUAL_C27_SAME_SOURCE_SIGNED_MU_ABEL_V1', 'ACTUAL_C27_SIGNED_JET_INTERFACE_V1',
            'ACTUAL_C27_STRUCTURED_TARGET_V1', 'ACTUAL_C27_FAR_COUPLING_GRAM_V1',
            'ACTUAL_C27_FINITE_DUAL_DEFECT_V1', 'ACTUAL_C27_FULL_GAMMA_PRIMAL_TAIL_V1',
            'ACTUAL_C27_SIGNED_HANKEL_TAIL_ENERGY_V1']:
        raise ValueError('Different original packet roles.')
    if (any(p['grid_bits'] != BITS for p in packets.values())
            or jet['original_metadata_grid_bits'] != 512 or far['bound_grid_bits'] != 628
            or defect['metadata'] != primal['metadata']
            or defect['metadata']['target_sha256'] != PINS['structured-joint-2048-target.json.zlib']
            or defect['metadata']['input_sha256'] != target['metadata']['input_sha256']
            or jet['input_sha256'] != target['metadata']['input_sha256']
            or (target['metadata']['N'], target['metadata']['M'], target['metadata']['last_mode'])
                != (FIRST-1, M, N)
            or target['metadata']['basis'] != defect['metadata']['basis']
            or target['metadata']['basis'] !=
                '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'
            or (far['c'], far['N'], far['M'], far['jet_order'], far['moment_diagonals'])
                != (27, FIRST-1, M, R, R)
            or (jet['physical_last_mode'], jet['first_far_mode'], jet['jet_order']) != (N, M+1, R)
            or jet['feature_order'] != 'A0..A55,B0..B55'
            or any(mu[k] != jet[k] or mu[k] != hankel[k] for k in
                   ('source_integer_sha256', 'source_bits', 'source_row'))
            or (mu['source_bits'], mu['source_row']) != (384, 0)):
        raise ValueError('Different physical support, source, phase, metadata grid or jet.')
    q, p = matrix_in(defect['D_inverse_difference'], BITS), matrix_in(defect['F_inverse_z'], BITS)
    h, kernel = matrix_in(jet['short_jet_coefficients'], BITS), matrix_in(target['kernel'], BITS)
    if ((q.nrows(), q.ncols(), p.nrows(), p.ncols()) != (N-FIRST+1, 1, N-FIRST+1, 1)
            or (h.nrows(), h.ncols(), kernel.nrows(), kernel.ncols()) != (2*R, 1, 2*R, 2*R)):
        raise ValueError('Incomplete original inverse profiles, feature direction or saved kernel.')
    scalar, rows, _, _ = inputs()  # Only read and merge the existing scalar archives.
    widths = [int(rows[n]['s'][1])-int(rows[n]['s'][0]) for n in range(N+1)]
    if scalar['grid_bits'] != BITS or min(widths) < 0 or max(widths) > 2:
        raise ValueError('The physical midpoint-radius premise fails.')
    midpoints = [arb(fmpq(int(rows[n]['s'][0])+int(rows[n]['s'][1]), 1 << (BITS+1)))
                 for n in range(FIRST, N+1)]
    left = left_features(list((q-p).entries()), midpoints, FIRST, M, R)
    rational = lambda packet, key: fmpq(*map(int, packet['exact_rationals'][key]))
    kappa, clin = rational(jet, 'kappa'), rational(jet, 'frozen_linear_envelope_upper')
    weight = rational(mu, 'original_far_weight')
    if (kappa != fmpq(33, 32) or kappa != rational(mu, 'kappa')
            or clin != rational(mu, 'frozen_linear_envelope_upper')
            or weight != fmpq(int(far['D_numerator']), 1 << far['D_grid_bits'])):
        raise ValueError('Different exact weight, legal coefficient or linear envelope.')
    symbol = interval_ball(primal['bounds']['actual_symbol_budget'], BITS)
    metadata_upper = lambda key: fmpq(int(target['metadata'][key][1]), 1 << 512)
    hankel_upper = lambda key: interval_ball(hankel['scalar_bounds'][key], BITS).upper()
    e0 = fmpq(int(far['zero_feature_pair'][1])-int(far['zero_feature_pair'][0]), 1 << (BITS+1))
    hstar = fmpq(int(far['upper_bounds']['moment_gram_norm']), 1 << far['bound_grid_bits'])
    exact, scalars, r2 = price_correction(symbol.lower(), symbol.upper(), weight,
        interval_ball(primal['bounds']['actual_mixed_tail_norm_upper'], BITS).upper(),
        metadata_upper('trial_operator_norm_sq_upper'), rational(jet, 'lift_norm_l1_upper'),
        metadata_upper('far_matrix_error'), hstar, e0, fmpq(1, 1 << BITS),
        hankel_upper('smooth_remainder_weighted_norm_upper'),
        rational(jet, 'selected_direction_norm_upper'),
        hankel_upper('principal_weighted_norm_upper'), kappa)
    nu0 = arb(kappa)*(left.transpose()*kernel*h)[0, 0]
    mu_interval = interval_ball(mu['scalar_bounds']['actual_proxy_functional'], BITS)
    allowance = exact['complete_theta_error_upper']
    theta = mu_interval+nu0+arb(0, allowance)
    delta = nu0+arb(0, allowance)
    center = arb(clin)+2*theta
    scalars.update(mixed_signed_kernel_cross=nu0, cached_actual_proxy_functional=mu_interval,
                   actual_minus_proxy_correction=delta, actual_signed_theta=theta,
                   center_upper_expression_before_good_metric=center)
    require_finite(list(scalars.values())+list(left.entries()))
    exports = {k: envelope(v, BITS) for k, v in scalars.items()}
    features = matrix_out(left, BITS)
    if any(v is None or int(v[0]) > int(v[1]) for v in exports.values()):
        raise ArithmeticError('Invalid signed-correction export; stop without retry.')
    exact.update(kappa=kappa, original_far_weight=weight, frozen_linear_envelope_upper=clin,
                 original_far_matrix_error=metadata_upper('far_matrix_error'),
                 physical_operator_norm_sq_upper=metadata_upper('trial_operator_norm_sq_upper'),
                 lift_norm_upper=rational(jet, 'lift_norm_l1_upper'), moment_gram_norm_upper=hstar,
                 zero_feature_radius=e0, scalar_midpoint_radius_upper=fmpq(1, 1 << BITS))
    packet = dict(schema='ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: digest for p, digest in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: digest for p, digest in guards.items() if p not in invoked},
        source_integer_sha256=mu['source_integer_sha256'], source_bits=384, source_row=0,
        first_correction_mode=FIRST, last_correction_mode=N, first_far_mode=M+1,
        normalization_scale=M, jet_order=R, feature_order=jet['feature_order'],
        grid_bits=BITS, working_bits=1024, original_metadata_grid_bits=512,
        moment_norm_grid_bits=far['bound_grid_bits'], fine_grid_bits=1024,
        fine_jet_remainder_squared=envelope(r2, 1024), scalar_bounds=exports,
        left_midpoint_features=features, physical_scalar_rows_checked=N+1,
        maximum_scalar_width_grid_units=max(widths),
        exact_rationals={k: [str(v.p), str(v.q)] for k, v in exact.items()},
        new_short_source_power_reductions=2, new_signed_kernel_crosses=1,
        new_far_features_or_moments=0, new_long_source_or_mu_evaluations=0,
        new_inverse_solve_trace_or_reference_calls=0, old_producers_replayed=0,
        parameter_source_precision_cutoff_scan=False, actual_signed_theta_interval_supplied=True,
        negative_upper_before_good_metric=bool(center.upper() < 0), good_block_metric_paid=False,
        actual_full_border_sign_claim=False, actual_rh_proof=False,
        verdict='conditional-same-source-signed-theta-interval',
        elapsed_seconds=time.monotonic()-started,
        conditions=mu['conditions']+[
            'The original physical Gram certificate applies to SAME Cphys on0..2048 with its original512-grid far_matrix_error, not the amplified complete allowance.',
            'q=D_inverse_difference and p=F_inverse_z; lambda=q-p on1281..2048, zero elsewhere; the original physical P0, saved K0 and cached h represent this common realization.',
            'Separate forward geometric/feature bounds use the true full moment norm, exact2/pi and actual symbol; the masked kernel is not assumed PSD.',
            'The signed correction replaces the old norm-only transfer allowance; the cached proxy interval already includes its source/symbol/Abel errors.'],
        scope='One fixed-source actual signed-theta enclosure under inherited actual-form/domain/source/lift and component certificates. No Lean, full-border, cofinal/retained positivity, Robin or RH conclusion.')
    verify_guards(guards)
    digest = write_packet(out, packet)
    print('saved', out.name, digest, flush=True)
    for key in ('mixed_signed_kernel_cross', 'sharp_forward_error_upper',
                'complete_theta_error_upper', 'actual_signed_theta',
                'center_upper_expression_before_good_metric'):
        print(key, scalars[key].str(24), flush=True)
    print('seconds', time.monotonic()-started, flush=True)


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
