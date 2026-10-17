"""Bound one actual far-row pairing using a zero-first-coordinate direction.

Reuse the saved signed response and source-energy upper. One new reference
action and two short-feature Gram contractions retain their common source.
No source, inverse, moment, sine or long-tail supplier is evaluated.
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
from full_gamma_coherent_split_capacity import freeze_upper
from full_gamma_correlated_sign import symmetric
from full_gamma_mixed_norm_reprice import mixed_norm
from full_gamma_negative_index_joint import freeze_lower
from full_gamma_shared_projection_energy import projected_energy
from full_gamma_signed_mu import require_finite, verify_guards, write_packet
from joint_gamma_boundary import matrix_in, matrix_out
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
PINS = {
    'full-gamma-directional-jet-norm.json.zlib':
        '9ccd6494c0520f9f53efe0b44bcae879566157c7e7ffb698165f9d02c258fd41',
    'full-gamma-relaxation-row-witness.json.zlib':
        '1879f4b0653e762ca0d7bb1e1933a59afc79fb10851d155f1c16c87055872992',
    'full-gamma-signed-jet-interface.json.zlib':
        '154949d8ab3e1f2875de531aa63a9c5a3c3d7d133640a5ce8b1d23c98d67565a',
    'full-gamma-negative-channel-capacity.json.zlib':
        '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
    'full-gamma-signed-correction.json.zlib':
        'd527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd',
    'joint-gamma-6765-rho32-v2-schur.json':
        '13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd',
    'full-gamma-directional-metric.json.zlib':
        'f6bea1b9b43c5ca8c85fc6cd7222038f2a40b81ed1ce2ef873af91ce1ca987a2',
    'full-gamma-projection-upper.json.zlib':
        'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3',
    'full-gamma-actual-channel-projection.json.zlib':
        '63e3c2e291940ca3549fa4a1303b55b4d92c08d594a50daf7ca2132087ad2745',
    'full-gamma-lift-projection.json.zlib':
        'a41813b27cd4567af326d42ec77e3b8c8942474b80bdf5833f7bc05aff52f873',
    'structured-joint-2048-target.json.zlib':
        'ea4a761f76b59a11b2d04faf80531901a5164f4145c2424a4e79889ae3f72e22',
}


def zero_direction(sigma, reference, rho, beta, alpha0, lam):
    """For P=T R T*, compute P*w*, w=rho-(rho[0]/lambda)*y0."""
    n = sigma.nrows()
    if (n < 2 or sigma.ncols() != n
            or (reference.nrows(), reference.ncols()) != (n-1, n-1)
            or (rho.nrows(), rho.ncols(), beta.nrows(), beta.ncols()) != (1, n, 1, n)):
        raise ValueError('Incompatible common finite-form dimensions.')
    require_finite(list(sigma.entries())+list(reference.entries())
                   +list(rho.entries())+list(beta.entries())+[arb(alpha0), arb(lam)])
    a = sigma[0, 0]
    if not a.upper() < 0 or lam.contains(0):
        raise ValueError('Missing negative pivot or nonzero reference coordinate.')
    y = arb(alpha0)*rho+beta
    if not y[0, 0].overlaps(lam):
        raise ValueError('Reference coordinate disagrees with the common finite row.')
    c = rho[0, 0]/lam
    # w[0]=0 is an exact common-parameter identity, not a rounded cancellation.
    wtail = arb_mat([[rho[0, j]-c*y[0, j] for j in range(1, n)]])
    transform = arb_mat([[-sigma[0, j+1]/a if i == 0 else arb(int(i == j+1))
                          for j in range(n-1)] for i in range(n)])
    direction = transform*(reference*wtail.transpose())
    require_finite([c]+list(wtail.entries())+list(direction.entries()))
    return c, wtail, direction


def source_pairing(features, direction, lift_diagonal, kernel, old_features,
                   eps_far, kphys, old_lift, old_norm2, source_energy, theta, kappa):
    """Apply existing physical Gram and common-projection bounds once."""
    n, rank = direction.nrows(), features.nrows()
    if (direction.ncols() != 1 or features.ncols() != n or n < 1 or rank < 1
            or len(lift_diagonal) != n or (kernel.nrows(), kernel.ncols()) != (rank, rank)
            or (old_features.nrows(), old_features.ncols()) != (rank, 1)):
        raise ValueError('Incompatible physical features, kernel or lift diagonal.')
    require_finite(list(features.entries())+list(direction.entries())
        +list(kernel.entries())+list(old_features.entries())+list(lift_diagonal)
        +[arb(x) for x in (eps_far, kphys, old_lift, old_norm2, source_energy, theta, kappa)])
    if (not kphys > 0 or not kappa > 0 or not old_norm2 > 0 or not source_energy > 0
            or min(eps_far, old_lift) < 0 or any(x.upper() < 0 for x in lift_diagonal)):
        raise ValueError('Missing physical-error, positive-norm or source-energy domain.')
    h = features*direction
    lift = freeze_upper(sum((abs(direction[j, 0]).upper()*lift_diagonal[j].upper().sqrt()
                            for j in range(n)), arb(0)).upper())
    physical_source_norm = freeze_upper(arb(kphys).sqrt()*arb(lift))
    quadratic, norm_error, norm2, _ = mixed_norm(h, kernel, eps_far, physical_source_norm)
    cross = (old_features.transpose()*kernel*h)[0, 0]
    cross_error = freeze_upper(arb(eps_far)*arb(kphys)*arb(old_lift)*arb(lift))
    cross_actual = cross+arb(0, arb(cross_error))
    # SAME u,b,d: <u,b>=theta/kappa, ||u||^2<=E, ||d||^2<=norm2.
    _, bounds = projected_energy(fmpq(0), source_energy, fmpq(0), norm2,
                                  old_norm2, theta/arb(kappa), cross_actual)
    lower = freeze_lower(arb(kappa)*arb(bounds['principal_mixed_cross_lower']))
    upper = freeze_upper(arb(kappa)*arb(bounds['principal_mixed_cross_upper']))
    scalars = dict(new_direction_gram_quadratic=quadratic,
        old_new_direction_kernel_cross=cross, old_new_direction_physical_cross=cross_actual)
    exact = dict(new_direction_lift_norm_upper=lift, new_direction_physical_norm2_upper=norm2,
        new_direction_norm_error_upper=norm_error, direction_cross_error_upper=cross_error,
        source_projection_pairing_lower=lower, source_projection_pairing_upper=upper)
    require_finite(list(scalars.values())+[arb(x) for x in exact.values()])
    return h, scalars, exact


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
            raise ValueError('Different accepted zero-channel input: '+name)
        guards[HERE/name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    latest, relaxed, jet, old, signed, core, cert, projection, actual, lift, target = [
        packets[name] for name in PINS]
    inherited = dict(PINS)
    for packet in (latest, relaxed, jet, signed):
        for field in ('invocation_source_sha256', 'inherited_dependency_sha256'):
            for name, digest in packet[field].items():
                if name in inherited and inherited[name] != digest:
                    raise ValueError('Conflicting inherited zero-channel dependency: '+name)
                inherited[name] = digest
    for name, digest in inherited.items():
        path = HERE/name
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked zero-channel dependency: '+name)
        guards[path] = digest
    verify_guards(guards)
    if ([p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_DIRECTIONAL_JET_NORM_V1', 'C27_SAVED_CONSTRAINT_RELAXATION_ROW_V1',
            'ACTUAL_C27_SIGNED_JET_INTERFACE_V1', 'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1',
            'ACTUAL_C27_SAME_SOURCE_SIGNED_CORRECTION_V1',
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1', 'ACTUAL_C27_FULL_GAMMA_DIRECTIONAL_METRIC_V1',
            'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1',
            'ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
            'ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1', 'ACTUAL_C27_STRUCTURED_TARGET_V1']
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any(p['source_integer_sha256'] != latest['source_integer_sha256']
                   for p in (jet, old, signed))
            or relaxed['reference_source_integer_sha256'] != latest['source_integer_sha256']
            or any((p['source_bits'], p['source_row']) != (384, 0) for p in (latest, jet, old, signed))
            or (relaxed['reference_source_bits'], relaxed['reference_source_row']) != (384, 0)
            or not relaxed['relaxed_negative_row_certified']
            or cert['packet_sha256']['core'] != PINS['joint-gamma-6765-rho32-v2-schur.json']
            or cert['positive_dimension'] != 17 or not cert['exact_reference_identities_checked']
            or cert['verdict'] != 'directional-center-enclosed'):
        raise ValueError('Different finite reference, scalar source or saved relaxation.')
    if (any(p['remaining_modes'] != old['remaining_modes']
            for p in (core, cert, projection, actual, lift))
            or len(old['remaining_modes']) != 18
            or any(p['input_sha256'] != old['input_sha256'] for p in (core, lift, jet))
            or target['metadata']['input_sha256'] != old['input_sha256']
            or lift['original_target_sha256'] != PINS['structured-joint-2048-target.json.zlib']
            or lift['original_metadata_grid_bits'] != 512
            or (lift['N'], lift['M'], lift['last_mode']) != (1280, 32768, 2048)
            or (jet['physical_last_mode'], jet['first_far_mode'], jet['jet_order'], jet['feature_order'])
                != (2048, 32769, 56, 'A0..A55,B0..B55')
            or lift['basis'] != target['metadata']['basis']
            or lift['basis'] != '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'):
        raise ValueError('Different physical lift, source, feature order or original grids.')
    rational = lambda p, key: fmpq(*map(int, p['exact_rationals'][key]))
    ball = lambda p, key: interval_ball(p['scalar_bounds'][key], BITS)
    energy, old_norm2 = [rational(latest, key) for key in
                         ('selected_far_energy_upper', 'selected_representer_energy_upper')]
    alpha0 = rational(old, 'fixed_alpha')
    kappa, weight = rational(latest, 'kappa'), rational(latest, 'original_far_weight')
    eps, kphys, old_lift = [rational(signed, key) for key in
        ('original_far_matrix_error', 'physical_operator_norm_sq_upper', 'lift_norm_upper')]
    if (kappa != fmpq(33, 32) or kappa != rational(signed, 'kappa')
            or weight != rational(signed, 'original_far_weight')
            or energy != rational(relaxed, 'source_energy_upper')
            or alpha0 != rational(relaxed, 'reference_alpha0')
            or old_lift != rational(jet, 'lift_norm_l1_upper')
            or eps != fmpq(int(target['metadata']['far_matrix_error'][1]), 1 << 512)
            or kphys != fmpq(int(target['metadata']['trial_operator_norm_sq_upper'][1]), 1 << 512)):
        raise ValueError('Different original physical Gram error or SAME source units.')
    # Metadata-only review stops before this first actual new direction conversion.
    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], BITS))
    reference = arb_mat([[arb(fmpq(x)) for x in row] for row in cert['reference_inverse']])
    rho = matrix_in([projection['matrices']['source_on_lift'][0]], BITS)
    beta = matrix_in(actual['actual_finite_row_on_lift'], BITS)
    c, wtail, direction = zero_direction(sigma, reference, rho, beta,
                                        alpha0, ball(old, 'finite_negative_channel'))
    features, kernel = matrix_in(lift['far_features_on_lift'], BITS), matrix_in(target['kernel'], BITS)
    old_features = matrix_in(jet['short_jet_coefficients'], BITS)
    diagonal = [interval_ball(lift['lift_gram'][j][j], BITS) for j in range(18)]
    if (features.nrows(), features.ncols(), kernel.nrows(), kernel.ncols()) != (112, 18, 112, 112):
        raise ValueError('Incomplete full56 physical feature/kernel input.')
    theta = ball(signed, 'actual_signed_theta')
    h, scalars, exact = source_pairing(features, direction, diagonal, kernel, old_features,
        eps, kphys, old_lift, old_norm2, energy, theta, kappa)
    rr, r0 = ball(old, 'steering_reference_rRr'), ball(old, 'selected_reference_form')
    p = ball(relaxed, 'reference_rho_positive_pairing')
    wenergy = rr-2*c*p+c*c*r0
    if wenergy.upper() < 0:
        raise ArithmeticError('Inconsistent positive zero-channel direction energy.')
    scalar_cap = freeze_upper((arb(rational(relaxed, 'relaxed_positive_energy_cap'))
                              *wenergy.upper()).sqrt())
    lo = max(-scalar_cap, exact['source_projection_pairing_lower'])
    hi = min(scalar_cap, exact['source_projection_pairing_upper'])
    if lo > hi:
        raise ArithmeticError('Source projection contradicts the retained scalar cap.')
    zeta = arb(arb((lo+hi)/2), arb((hi-lo)/2))
    # Q_R*w*=P*w* since w[0]=0; rho=w+c*y0 on this same finite core.
    derivative = zeta+c*theta
    s = rational(relaxed, 'fixed_relaxation_alpha')-alpha0
    relaxed_zeta = ball(relaxed, 'relaxed_row_coefficient')*(p-c*r0)-arb(s)*(rr-c*p)
    disjoint = bool(zeta.upper() < relaxed_zeta.lower() or relaxed_zeta.upper() < zeta.lower())
    exact.update(source_energy_upper=energy, old_representer_norm2_upper=old_norm2,
        scalar_only_zero_pairing_absolute_upper=scalar_cap,
        zero_pairing_lower=lo, zero_pairing_upper=hi, original_far_weight=weight, kappa=kappa,
        zero_pairing_absolute_upper=max(abs(lo), abs(hi)))
    scalars.update(coordinate_ratio=c, zero_direction_positive_energy=wenergy,
        actual_far_zero_pairing=zeta, actual_far_reference_rho_pairing=derivative,
        constructed_relaxation_zero_pairing=relaxed_zeta)
    scalars.update({key: arb(value) for key, value in exact.items()})
    require_finite(list(scalars.values()))
    exported = {key: envelope(value, BITS) for key, value in scalars.items()}
    arrays = dict(zero_direction_tail=matrix_out(wtail, BITS),
                  positive_reference_direction_on_lift=matrix_out(direction, BITS),
                  new_short_jet_coefficients=matrix_out(h, BITS))
    if (any(pair is None or int(pair[0]) > int(pair[1]) for pair in exported.values())
            or any(pair is None or int(pair[0]) > int(pair[1])
                   for rows in arrays.values() for row in rows for pair in row)):
        raise ArithmeticError('Invalid zero-channel interval export; stop without retry.')
    result = dict(schema='ACTUAL_C27_ZERO_CHANNEL_PAIRING_BOUND_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            path.name: digest for path, digest in guards.items() if path in invoked},
        inherited_dependency_sha256={path.name: digest for path, digest in guards.items() if path not in invoked},
        input_sha256=old['input_sha256'], remaining_modes=old['remaining_modes'],
        source_integer_sha256=latest['source_integer_sha256'], source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=1024, original_metadata_grid_bits=512,
        exact_rationals={key: [str(value.p), str(value.q)] for key, value in exact.items()},
        scalar_bounds=exported, **arrays,
        source_pairing_strictly_narrower_than_scalar_cap=bool(max(abs(lo), abs(hi)) < scalar_cap),
        constructed_relaxation_row_excluded_by_source_pairing=disjoint,
        actual_far_pairing_interval_supplied=True, actual_derivative_value_supplied=False,
        new_reference_vector_actions=1, new_cached_feature_actions=1,
        new_physical_feature_quadratics=1, new_physical_feature_crosses=1,
        old_direction_quadratic_replayed=False,
        new_inverse_solve_trace_source_sine_moment_or_tail_suppliers=0,
        old_producers_replayed=0, parameter_source_precision_cutoff_scan=False,
        actual_center_sign_claim=False, actual_full_border_sign_claim=False, actual_rh_proof=False,
        conditions=latest['conditions']+[
            'On the SAME finite core, c=rho[0]/lambda, y0=alpha0*rho+beta and w=rho-c*y0 have w[0]=0 exactly. Thus Q_R*w*=P*w* and f Q_R rho*=f P w*+c theta.',
            'The actual far row satisfies f v=kappa*<u_F,Bv> in the SAME original weighted Hilbert space. Therefore f P w*=kappa*<u_F,d>, d=B P w*, and theta=kappa*<u_F,b>, b=B Q_R y0*.',
            'The accepted original physical Gram operator error pays the new norm by eps_far*Kphys*Ld^2 and the old/new cross by eps_far*Kphys*Lb*Ld. The old b norm is the saved sharper cap; its quadratic is not replayed. Rounded or masked kernel centers need not be PSD.',
            'The existing common-representer Gram bound uses the SAME u_F,b,d, signed theta interval and source-energy UPPER; no signed pairing is inferred from a norm alone. The scalar-only P cap is intersected in the same direction.',
            'The constructed relaxation pairing retains its own saved tau and chosen theta, not the actual theta. Disjoint intervals would exclude that particular row under the added source-energy/Gram relation, not all negative rows or the actual paid expression.',
            'This supplies conditional intervals, not an exact derivative, actual paid comparison, physical steering legality, center/full-border/cofinal positivity, Robin or RH.'
        ], scope='One numerical application of existing linear, physical-Gram and common-projection bounds. No generic theorem, novelty or Lean claim.',
        elapsed_seconds=time.monotonic()-started)
    verify_guards(guards)
    digest = write_packet(out, result)
    print('saved', out.name, digest, flush=True)
    for key in ('actual_far_zero_pairing', 'constructed_relaxation_zero_pairing',
                'actual_far_reference_rho_pairing', 'scalar_only_zero_pairing_absolute_upper'):
        print(key, scalars[key].str(24), flush=True)
    print('narrower', result['source_pairing_strictly_narrower_than_scalar_cap'],
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
