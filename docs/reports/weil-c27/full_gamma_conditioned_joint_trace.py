"""Bound the complete saved joint-border budget with conditional signed traces.

Only saved-data regression, gradient and signed trace contractions run.
The finite reference diagonal and old coordinate-box bound are reused. No
source query, old scientific producer, inverse, solve or scan is replayed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import zlib

from flint import arb, arb_mat, ctx, fmpq
from diagonal_floor_witness import interval_ball
from full_gamma_correlated_sign import symmetric
from full_gamma_far_energy_budget import finite
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
PINS = {
    'full-gamma-conditioned-profile.json.zlib': '2e1c11ef3807445ce6b157d47232bc6086ac96f8a07786e4b03b7355e61c00ae',
    'joint-gamma-6765-rho32-v2-schur.json': '13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd',
    'full-gamma-directional-metric.json.zlib': 'f6bea1b9b43c5ca8c85fc6cd7222038f2a40b81ed1ce2ef873af91ce1ca987a2',
    'full-gamma-negative-channel-capacity.json.zlib': '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
    'full-gamma-actual-channel-projection.json.zlib': '63e3c2e291940ca3549fa4a1303b55b4d92c08d594a50daf7ca2132087ad2745',
    'full-gamma-far-energy-budget.json.zlib': '11ccabd6110436b2b29f1a910c6c36a5fd8479af19d183cf614826ed9ad746c6',
    'full-gamma-projection-upper.json.zlib': 'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3',
    'full-gamma-finite-good-block-metric.json.zlib': '9ba53b70a0f1cad25b9e92894eff93f1340dc31d1a9bf945c560b6dc6940d072',
    'full-gamma-joint-far-row.json.zlib': '00b9e4bdd8a99da4d6fdafc696d313bd43e4060e1ba7d8cec316f7cc5bcb7707',
}
PROFILE, CORE, CERT, CAPACITY, FINITE, BUDGET, PROJECTION, METRIC, FAR = PINS


def conditioned_energy(sigma, reference, beta, regression, measured,
                       conditional, remaining, eta, finite00, old_upper):
    """PSD trace monotonicity and Minkowski, for one common realization.

    With p<0, T=[-Sigma[0,1:]/p;I], R>0 and R/(1-eta) >= Pplus^-1,
    r=y-mB obeys r^T*r <= remaining*S_F. Thus its reference energy is
    at most remaining*tr(R*T^T*S_F*T). The finite00 is an already paid
    beta_J*T*R*T^T*beta_J^T enclosure, not a new reference action.
    """
    n = sigma.nrows()
    if (n < 2 or sigma.ncols() != n
            or (reference.nrows(), reference.ncols()) != (n-1, n-1)
            or (beta.nrows(), beta.ncols()) != (1, n)
            or (regression.nrows(), regression.ncols()) != (2, n)
            or (measured.nrows(), measured.ncols()) != (1, 2)
            or (conditional.nrows(), conditional.ncols()) != (n, n)):
        raise ValueError('Incomplete common-source conditional core domains.')
    finite(list(sigma.entries())+list(reference.entries())+list(beta.entries())
           +list(regression.entries())+list(measured.entries())+list(conditional.entries())
           +[arb(remaining), arb(eta), finite00, arb(old_upper)])
    pivot = sigma[0, 0]
    if (not pivot < 0 or not 0 <= eta < 1 or not remaining >= 0
            or not old_upper >= 0 or finite00.upper() < 0):
        raise ArithmeticError('Missing pivot, inherited reference or energy domain.')
    transport = arb_mat([[-sigma[0, j+1]/pivot if i == 0 else arb(int(i == j+1))
                         for j in range(n-1)] for i in range(n)])
    bt, rt = beta*transport, regression*transport
    actions = reference*rt.transpose()
    mixed, reg_gram = bt*actions, symmetric(rt*actions)
    mean_gram = arb_mat(3, 3)
    mean_gram[0, 0] = finite00
    for i in range(2):
        mean_gram[0, i+1] = mean_gram[i+1, 0] = mixed[0, i]
        for j in range(2):
            mean_gram[i+1, j+1] = reg_gram[i, j]
    coefficients = arb_mat([[1, measured[0, 0], measured[0, 1]]])
    mean_energy = (coefficients*mean_gram*coefficients.transpose())[0, 0]
    conditional_positive = symmetric(transport.transpose()*conditional*transport)
    weighted = reference*conditional_positive
    trace = sum((weighted[i, i] for i in range(n-1)), arb(0))
    if trace.upper() < 0 or mean_energy.upper() < 0:
        raise ArithmeticError('Inherited PSD reference geometry has no nonnegative enclosure.')
    residual_upper = arb(remaining)*max(arb(0), trace.upper())
    mean_upper = max(arb(0), mean_energy.upper())
    reference_upper = (mean_upper.sqrt()+residual_upper.sqrt())**2
    actual_upper = reference_upper/(1-arb(eta))
    scalars = dict(reused_finite_reference_diagonal=finite00,
        mean_reference_energy=mean_energy, conditional_reference_trace=trace,
        conditional_residual_reference_energy_upper=residual_upper,
        complete_reference_energy_upper=reference_upper,
        actual_positive_core_energy_upper=actual_upper)
    finite(list(mean_gram.entries())+list(conditional_positive.entries())
           +list(actions.entries())+list(scalars.values()))
    new_upper = fmpq(int(envelope(actual_upper, BITS)[1]), 1 << BITS)
    exact = dict(remaining_energy_upper=remaining, neumann_eta_upper=eta,
        original_grid_positive_core_energy_upper=old_upper,
        new_grid_positive_core_energy_upper=new_upper)
    if old_upper > 0:
        exact['new_to_old_upper_ratio'] = new_upper/old_upper
    result = dict(regression_on_positive_core=matrix_out(rt, BITS),
        reference_regression_actions=matrix_out(actions, BITS),
        mean_reference_gram=matrix_out(mean_gram, BITS),
        conditional_positive_gram=matrix_out(conditional_positive, BITS),
        scalar_bounds={k: envelope(v, BITS) for k, v in scalars.items()},
        exact_rationals={k: [str(v.p), str(v.q)] for k, v in exact.items()},
        verdict=('strictly_tighter_positive_core_energy_upper' if new_upper < old_upper
                 else 'no_positive_core_energy_upper_improvement'))
    return result, transport, actions, mean_gram, reference_upper


def joint_trace_budget(sigma, reference, beta, regression, measured, conditional,
                       remaining, eta, finite00, old_upper, rho, gram, metric_schur,
                       border_a, charge, saved_rr, saved_rb, raw_far_pairs):
    """A sufficient upper for the FULL zcrit-QJ, not separate extreme terms.

    The nominal inverse uses the paid reference R. Its conditional curvature
    is bounded by Wolkowicz--Styan Theorem2.1, Eq(2.3), pp474-476:
    https://www.math.uwaterloo.ca/~hwolkowi/henry/reports/bndseigs80.pdf
    PSD norm estimates
    price the entire actual-minus-nominal expression, including the saved
    Neumann error and the actual border denominator. No S inverse is used.
    """
    result, t, actions, mean_gram, reference_upper = conditioned_energy(
        sigma, reference, beta, regression, measured, conditional,
        remaining, eta, finite00, old_upper)
    n = sigma.nrows()
    if ((rho.nrows(), rho.ncols()) != (1, n)
            or (gram.nrows(), gram.ncols()) != (2, 2)
            or len(raw_far_pairs) != n
            or any(len(p) != 2 or int(p[0]) > int(p[1]) for p in raw_far_pairs)):
        raise ValueError('Incomplete joint border, metric or original row box.')
    finite(list(rho.entries())+list(gram.entries())
           +[metric_schur, border_a, arb(charge), saved_rr, saved_rb])
    if (not gram[0, 0] > 0 or not metric_schur > 0 or not border_a > 0
            or not charge >= 0 or saved_rr.upper() < 0):
        raise ArithmeticError('Missing positive observed metric or border denominator.')
    a_pair = envelope(border_a, BITS)
    a0 = fmpq(int(a_pair[0])+int(a_pair[1]), 1 << (BITS+1))
    p = sigma[0, 0]
    wmetric = symmetric(t*reference*t.transpose())
    ghat = wmetric + arb_mat([[1/p if i == j == 0 else 0
                              for j in range(n)] for i in range(n)])
    rhat = ghat*rho.transpose()
    h = arb_mat(1, 3)
    h[0, 0] = gram[0, 1]+rho[0, 0]*beta[0, 0]/p+saved_rb
    rho_regression = (rho*t)*actions
    for j in range(2):
        h[0, j+1] = (int(j == 0)+rho[0, 0]*regression[j, 0]/p
                     +rho_regression[0, j])
    v = arb_mat([[beta[0, i], regression[0, i], regression[1, i]] for i in range(n)])
    first = arb_mat([[v[0, j] for j in range(3)]])
    mean_signed = mean_gram+first.transpose()*first/p
    c, s, a = gram[0, 1], metric_schur, gram[0, 0]
    qj_inverse = arb_mat([[1/a+c*c/(a*a*s), -c/(a*s)], [-c/(a*s), 1/s]])
    cq = h.transpose()*h/arb(a0)-mean_signed
    cq[0, 0] -= arb(charge)+gram[1, 1]
    cq[0, 2] -= 1
    cq[2, 0] -= 1
    for i in range(2):
        for j in range(2):
            cq[i+1, j+1] -= qj_inverse[i, j]
    cq = symmetric(cq)
    z = ghat*(rho.transpose()*h/arb(a0)-v)
    cg = symmetric(z.transpose()*conditional*z)
    hb = symmetric(rhat*rhat.transpose()/arb(a0)-ghat)
    sh = conditional*hb
    t1 = sum((sh[i, i] for i in range(n)), arb(0))
    sh2 = sh*sh
    t2 = sum((sh2[i, i] for i in range(n)), arb(0))
    distance = max(arb(0), abs(t1).lower())
    discriminant = (t2.upper()/n-distance**2/(n*n)).upper()
    if t2.upper() < 0 or discriminant < 0:
        raise ArithmeticError('Signed trace enclosures contradict the inherited PSD geometry.')
    curvature = max(arb(0), (t1.upper()/n+((n-1)*discriminant).sqrt()).upper())
    theta = arb_mat([[1], [measured[0, 0]], [measured[0, 1]]])
    qc = (theta.transpose()*cq*theta)[0, 0]
    gradient_energy = (theta.transpose()*cg*theta)[0, 0]
    qnorm = (rhat.transpose()*conditional*rhat)[0, 0]
    if gradient_energy.upper() < 0 or qnorm.upper() < 0:
        raise ArithmeticError('Conditional gradient Gram has no PSD realization.')
    ellipsoid_linear = (arb(remaining)*max(arb(0), gradient_energy.upper())).sqrt().upper()
    gradient = z*theta
    mu = regression.transpose()*measured.transpose()
    midpoint = arb_mat([[arb(fmpq(int(lo)+int(hi), 1 << (BITS+1)))]
                       for lo, hi in raw_far_pairs])
    radii = [fmpq(int(hi)-int(lo), 1 << (BITS+1)) for lo, hi in raw_far_pairs]
    box_linear = ((gradient.transpose()*(midpoint-mu))[0, 0]
                  +sum((arb(radii[i])*abs(gradient[i, 0]) for i in range(n)), arb(0)))
    linear_upper = min(ellipsoid_linear, box_linear.upper())
    affine_upper = qc.upper()+2*linear_upper
    nominal_upper = affine_upper+arb(remaining)*curvature
    relative = arb(eta)/(1-arb(eta))
    ref_cap, rr_cap = max(arb(0), reference_upper.upper()), max(arb(0), saved_rr.upper())
    delta_o = relative*(rr_cap*ref_cap).sqrt()
    ohat_upper = (abs((h*theta)[0, 0]).upper()
                  +(arb(remaining)*max(arb(0), qnorm.upper())).sqrt().upper())
    denominator_error = abs(1/border_a-1/arb(a0)).upper()
    perturbation = ((2*ohat_upper*delta_o+delta_o**2)/border_a.lower()
                    +ohat_upper**2*denominator_error+relative*ref_cap)
    upper = nominal_upper+perturbation
    values = dict(nominal_signed_trace=t1, nominal_signed_square_trace=t2,
        nominal_curvature_upper=curvature, nominal_center_quadratic=qc,
        nominal_gradient_energy=gradient_energy, nominal_source_direction_energy=qnorm,
        ellipsoid_linear_upper=arb(ellipsoid_linear), original_box_linear_upper=box_linear,
        selected_linear_upper=arb(linear_upper), affine_upper=affine_upper,
        nominal_budget_upper=nominal_upper, reference_error_factor=relative,
        nominal_offdiagonal_magnitude_upper=ohat_upper, offdiagonal_reference_error_upper=delta_o,
        inverse_denominator_error_upper=arb(denominator_error),
        complete_reference_perturbation_upper=perturbation, complete_joint_budget_upper=upper)
    finite(list(ghat.entries())+list(cq.entries())+list(cg.entries())
           +list(z.entries())+list(values.values()))
    result['positive_energy_verdict'] = result.pop('verdict')
    result['scalar_bounds'].update({k: envelope(v, BITS) for k, v in values.items()})
    result.update(nominal_center_quadratic_matrix=matrix_out(cq, BITS),
        nominal_gradient=matrix_out(z, BITS), nominal_gradient_gram=matrix_out(cg, BITS),
        nominal_border_quadratic=matrix_out(hb, BITS))
    upper_q = fmpq(int(envelope(upper, BITS)[1]), 1 << BITS)
    result['exact_rationals'].update(nominal_border_denominator=[str(a0.p), str(a0.q)],
        complete_joint_budget_upper=[str(upper_q.p), str(upper_q.q)])
    result['verdict'] = ('fixed_optimistic_comparison_excluded_by_joint_trace'
                         if upper_q < 0 else 'inconclusive_conditional_joint_trace_budget')
    if upper_q < 0:
        determinant_lower = fmpq(int(a_pair[0]), 1 << BITS)*(-upper_q)
        result['exact_rationals']['fixed_comparison_determinant_positive_lower'] = [
            str(determinant_lower.p), str(determinant_lower.q)]
    return result


def consume(out):
    started = time.perf_counter()
    packets = {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Input hash mismatch: '+name)
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    profile, core, cert, capacity, actual, budget, projection, metric, far = packets.values()
    if [p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_CONDITIONED_TWENTY_PROFILE_GRAM_V1',
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1',
            'ACTUAL_C27_FULL_GAMMA_DIRECTIONAL_METRIC_V1',
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1',
            'ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
            'ACTUAL_C27_FAR_INVERSE_ENERGY_BUDGET_V1',
            'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1',
            'ACTUAL_C27_FINITE_GOOD_BLOCK_JOINT_METRIC_V1',
            'ACTUAL_C27_JOINT_FAR_ROW_V1']:
        raise ValueError('Different conditional/core/reference/source packet roles.')
    if (not profile['conditioning_certified'] or not cert['exact_reference_identities_checked']
            or cert['positive_dimension'] != 17
            or any(p['remaining_modes'] != core['remaining_modes']
                   for p in (profile, cert, capacity, actual, budget, projection, far))
            or len(core['remaining_modes']) != 18
            or any(p['grid_bits'] != BITS for p in packets.values())
            or profile['original_metadata_grid_bits'] != 512
            or any((p['source_integer_sha256'], p['source_bits'], p['source_row'])
                   != (SOURCE_ID, 384, 0) for p in (profile, capacity, budget, metric, far))):
        raise ValueError('Different source, original lift, grids or inherited conditioning.')
    for p in (capacity, budget):
        if any(p['packet_sha256'][n] != PINS[n] for n in (CORE, CERT, FINITE)):
            raise ValueError('Different inherited reference and finite row.')
    if (budget['packet_sha256'][CAPACITY] != PINS[CAPACITY]
            or cert['packet_sha256']['core'] != PINS[CORE]
            or profile['packet_sha256']['full-gamma-lift-projection.json.zlib']
            != capacity['packet_sha256']['full-gamma-lift-projection.json.zlib']
            or actual['input_sha256'] != core['input_sha256']
            or cert['input_sha256'] != core['input_sha256']):
        raise ValueError('Different same-target/source/lift normalization.')
    if (any(budget['packet_sha256'][n] != PINS[n] for n in (PROJECTION, METRIC, FAR))
            or any(profile['packet_sha256'][n] != PINS[n] for n in (METRIC, FAR))
            or capacity['packet_sha256'][PROJECTION] != PINS[PROJECTION]
            or metric['basis'] != profile['basis']
            or metric['positive_modes'] != profile['positive_modes']):
        raise ValueError('Different observed inverse metric, original row or projection.')
    helpers = profile['inherited_dependency_sha256']
    for name, digest in helpers.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Inherited helper hash mismatch: '+name)
    alpha = fmpq(*map(int, capacity['exact_rationals']['fixed_alpha']))
    eta = fmpq(int(cert['scalar_bounds']['neumann_eta_upper'][1]), 1 << BITS)
    if eta != fmpq(*map(int, budget['exact_rationals']['neumann_eta_upper'])):
        raise ValueError('Different inherited Neumann allowance.')
    bounds = capacity['scalar_bounds']
    finite00 = (interval_ball(bounds['selected_reference_form'], BITS)
                -2*arb(alpha)*interval_ball(bounds['steering_reference_rRb'], BITS)
                -arb(alpha)**2*interval_ball(bounds['steering_reference_rRr'], BITS))
    remaining = fmpq(*map(int, profile['exact_rationals']['remaining_energy_upper']))
    old = budget['scalar_bounds']
    old_upper = fmpq(int(old['reference_complete_energy'][1])
                     +int(old['reference_bb_error'][1]), 1 << BITS)
    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], BITS))
    reference = arb_mat([[arb(fmpq(v)) for v in row] for row in cert['reference_inverse']])
    result = joint_trace_budget(sigma, reference,
        matrix_in(actual['actual_finite_row_on_lift'], BITS),
        matrix_in(profile['regression'], BITS), matrix_in(profile['measured_pair'], BITS),
        symmetric(matrix_in(profile['conditional_gram'], BITS)), remaining, eta,
        finite00, old_upper,
        matrix_in([projection['matrices']['source_on_lift'][0]], BITS),
        symmetric(matrix_in(metric['joint_inverse_metric'], BITS)),
        interval_ball(budget['scalar_bounds']['metric_schur'], BITS),
        interval_ball(budget['scalar_bounds']['border_a'], BITS),
        fmpq(*map(int, budget['exact_rationals']['fixed_charge'])),
        interval_ball(bounds['steering_reference_rRr'], BITS),
        interval_ball(bounds['steering_reference_rRb'], BITS), far['actual_far_row'][0])
    result.update(schema='ACTUAL_C27_CONDITIONED_JOINT_TRACE_BUDGET_V1',
        packet_sha256=PINS, inherited_dependency_sha256=helpers,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_integer_sha256=SOURCE_ID, source_bits=384, source_row=0,
        remaining_modes=core['remaining_modes'], grid_bits=BITS, working_bits=ctx.prec,
        original_metadata_grid_bits=profile['original_metadata_grid_bits'],
        exact_far_allowance=profile['exact_far_allowance'],
        source_energy_upper=profile['source_energy_upper'], kappa=profile['kappa'],
        comparison_t_plus=budget['exact_rationals']['comparison_t_plus'],
        fixed_charge=budget['exact_rationals']['fixed_charge'],
        comparison_coefficient_eligibility_claim=False,
        new_regression_reference_vector_actions=2, new_mixed_regression_gram_entries=5,
        new_conditional_positive_gram_transport=1, new_conditional_reference_trace=1,
        new_nominal_core_metric_contraction=1, new_signed_trace_contractions=2,
        new_gradient_gram_dimension=3,
        reused_finite_reference_diagonal=True, reused_original_grid_baseline=True,
        new_source_far_target_inverse_solve_queries=0,
        old_scientific_producers_or_fixtures_replayed=0,
        source_direction_parameter_precision_cutoff_scan=False,
        full_entry_box_asserted_PSD=False, unseen_far_inverse_energy_upper_supplied=False,
        actual_center_sign_claim=False, actual_full_border_sign_claim=False,
        actual_rh_proof=False, elapsed_seconds=time.perf_counter()-started,
        inherited_profile_conditions=profile['conditions'],
        conditions=[
            'All inherited common-form/domain, physical-majorant, source-energy, source/lift and normalization premises in inherited_profile_conditions remain required. Parent-result limitations there describe the parent packet.',
            'The SAME finite beta_J, measured m, actual far row y=fV and exact lift are used. No endpoints are asserted jointly attainable.',
            'The inherited actual physical majorant gives exact S_F PSD, not signed-kernel or full interval-entry-box PSD. Its original epsilon/512-metadata allowance is unchanged.',
            'The accepted exact R is positive and R/(1+eta)<=Pplus^-1<=R/(1-eta), with the original eta<1. No reference=actual-inverse identity is asserted.',
            'The residual satisfies r^T*r <= remaining*S_F, including the saved source energy and raw kappa exactly once. PSD trace monotonicity bounds its R-energy after transport.',
            'The finite reference diagonal is reused through the saved identity gRg-2*alpha*rRb-alpha^2*rRr; only the new two regression rows receive reference actions.',
            'Minkowski combines a uniform mean energy upper and residual energy upper for one realization. The original paid row-box baseline is read from raw640-grid bounds, without round-trip widening.',
            'The complete nominal zcrit-QJ is assembled before interval evaluation. Its signed curvature uses tr(S*H) and tr(S*H*S*H), applying the existing real-eigenvalue trace bound in dimension18. No positive-kernel or full entry-box PSD hypothesis is introduced.',
            'The nominal gradient linear term is bounded by the smaller of a conditional-ellipsoid support and the original raw640-grid coordinate-box support; no marginal shrink or simultaneous endpoint attainment is assumed.',
            'The entire actual-minus-nominal budget is bounded using eta/(1-eta), the NEW complete conditional reference-energy upper, saved rho-reference energy, and the saved actual positive border-denominator enclosure. No reference or midpoint metric is substituted for the actual inverse without this payment.',
            'A strictly negative upper for zcrit-QJ gives det(N2)>0 using the EXISTING LOWER bound z>=QJ; it excludes only the fixed optimistic comparison represented here. It does not assert parameter eligibility, target transfer or a z upper. An inconclusive upper is not a counterexample or a proof of impossibility.',
            'Actual-center/full-border/cofinal positivity, Robin and RH remain unresolved. Existing trace, projection, norm and Neumann methods are applied; no generic-method novelty or Lean result is claimed.'])
    out.parent.mkdir(parents=True, exist_ok=True)
    write_packet(out, result)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output; there is no source or precision retry.')
    ctx.prec = 1024
    consume(args.out)


if __name__ == '__main__':
    main()
