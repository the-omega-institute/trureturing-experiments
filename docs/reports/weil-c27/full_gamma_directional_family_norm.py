"""Use source-response geometry for the same exact negative family.

Reuse the normalized primal margin, coefficient ellipse, conditional source
geometry, relative inverse certificate and original allowance. Only the
five-response Euclidean Gram and conditional Euclidean trace are new.
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
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
FAMILY = 'full-gamma-coherent-family-norm.json.zlib'
FAMILY_SHA = '6d81b66716ed241500982bdc3e46cf5b096a1e4922759b1ff26b255f763b2018'
PROFILE = 'full-gamma-conditioned-profile.json.zlib'
CORE = 'joint-gamma-6765-rho32-v2-schur.json'
CERT = 'full-gamma-directional-metric.json.zlib'
METRIC = 'full-gamma-finite-good-block-metric.json.zlib'
BUDGET = 'full-gamma-far-energy-budget.json.zlib'
TRACE = 'full-gamma-conditioned-joint-trace.json.zlib'
ACTUAL = 'full-gamma-actual-channel-projection.json.zlib'
PROJECTION = 'full-gamma-projection-upper.json.zlib'
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'


def finite(values):
    if any(not arb(v).is_finite() for v in values):
        raise ArithmeticError('Nonfinite directional family data.')


def directional_geometry(sigma, reference, euclidean, gram, schur,
                         rho, beta, regression_actions, conditional):
    """Exact enclosed responses; saved regression actions are not replayed."""
    n = sigma.nrows()
    dims = [(sigma, n, n), (reference, n-1, n-1),
            (euclidean, n+2, n+2), (gram, 2, 2), (rho, 1, n),
            (beta, 1, n), (regression_actions, n-1, 2), (conditional, n, n)]
    if n < 2 or any((x.nrows(), x.ncols()) != (rows, cols) for x, rows, cols in dims):
        raise ValueError('Different SAME-frame directional geometry dimensions.')
    finite([v for x, _, _ in dims for v in x.entries()]+[schur])
    pivot, a, c = sigma[0, 0], gram[0, 0], gram[0, 1]
    if not pivot < 0 or not a > 0 or not schur > 0:
        raise ArithmeticError('Missing actual negative pivot or positive good metric.')
    inverse_good = arb_mat([[1/a+c*c/(a*a*schur), -c/(a*schur)],
                           [-c/(a*schur), 1/schur]])
    transport = arb_mat([[-sigma[0, j+1]/pivot if i == 0 else arb(int(i == j+1))
                          for j in range(n-1)] for i in range(n)])
    # Two previously unsaved response products; no inverse, solve or old consume.
    two_rows = arb_mat([[rho[0, j] for j in range(n)],
                       [beta[0, j] for j in range(n)]])
    new_actions = reference*(two_rows*transport).transpose()
    actions = arb_mat(n-1, 4)
    for i in range(n-1):
        for j in range(2):
            actions[i, j] = new_actions[i, j]
            actions[i, j+2] = regression_actions[i, j]
    lifted = transport*actions
    responses = arb_mat(n+2, 5)
    responses[2, 0] = 1
    responses[0, 1] = responses[1, 2] = 1
    for i in range(2):
        for j in range(2):
            responses[i, j+3] = inverse_good[i, j]
    for i in range(n):
        for j in range(4):
            responses[i+2, j+1] = lifted[i, j]
    weighted = euclidean*responses
    response_gram = arb_mat(5, 5)
    response_gram[0, 0] = euclidean[2, 2]  # Saved negative-column norm, not a new dot.
    for i in range(5):
        for j in range(i, 5):
            if i == j == 0:
                continue
            value = sum((responses[k, i]*weighted[k, j] for k in range(n+2)), arb(0))
            response_gram[i, j] = response_gram[j, i] = value
    vv = arb_mat([[euclidean[i+2, j+2] for j in range(n)] for i in range(n)])
    response_kernel = symmetric(transport*reference*transport.transpose())
    product = vv*response_kernel*conditional*response_kernel
    residual_trace = sum((product[i, i] for i in range(n)), arb(0))
    if residual_trace.upper() < 0:
        raise ArithmeticError('Conditional PSD response trace has a negative upper.')
    finite(list(response_gram.entries())+[residual_trace])
    return response_gram, residual_trace


def directional_budget(response_gram, measured, negative_norm, margin,
                       energy, t_plus, fee, remaining, eta, core_trace,
                       residual_trace, delta):
    """The SAME exact p satisfies p*Lambda^-1*p <= energy, even with q."""
    if ((response_gram.nrows(), response_gram.ncols()) != (5, 5)
            or (measured.nrows(), measured.ncols()) != (1, 2)):
        raise ValueError('Incomplete five-response family or measured pair.')
    scalars = (negative_norm, margin, energy, t_plus, fee, remaining,
               eta, core_trace, residual_trace, delta)
    finite(list(response_gram.entries())+list(measured.entries())+list(scalars))
    if (not margin > 0 or not t_plus > 0 or not fee > 0 or not 0 <= eta < 1
            or any(not x >= 0 for x in (negative_norm, energy, remaining,
                                       core_trace, residual_trace, delta))):
        raise ArithmeticError('Missing exact-family ellipse, trace or allowance domain.')
    coefficients = arb_mat([[1], [measured[0, 0]], [measured[0, 1]]])
    mean_gram = arb_mat([[response_gram[i+2, j+2] for j in range(3)] for i in range(3)])
    mean_norm = (coefficients.transpose()*mean_gram*coefficients)[0, 0]
    cross = sum((response_gram[1, j+2]*coefficients[j, 0] for j in range(3)), arb(0))
    negative_cross = sum((response_gram[0, j+2]*coefficients[j, 0] for j in range(3)), arb(0))
    if mean_norm.upper() < 0 or response_gram[1, 1].upper() < 0:
        raise ArithmeticError('Response PSD geometry contradicts its norm enclosure.')
    aa = (arb(t_plus)*max(arb(0), response_gram[1, 1].upper())).upper()
    dd = (max(arb(0), mean_norm.upper())/arb(fee)).upper()
    cc = ((arb(t_plus)/arb(fee)).sqrt()*abs(cross)).upper()
    # Rayleigh comparison by absolute offdiagonal; this is an UPPER eigenvalue.
    top = ((aa+dd+((aa-dd)**2+4*cc**2).sqrt())/2).upper()
    cross_energy = (arb(t_plus)*response_gram[0, 1]**2
                    +negative_cross**2/arb(fee)).upper()
    if cross_energy < 0:
        raise ArithmeticError('Negative source-cross support upper.')
    nominal = (arb(negative_norm)+arb(energy)*top
               +2*(arb(energy)*cross_energy).sqrt()).upper()
    residual = (arb(energy)/arb(fee)*arb(remaining)*arb(residual_trace)).sqrt().upper()
    reference_error = (arb(eta)/(1-arb(eta))
        *(arb(core_trace)*(1+arb(eta))*arb(energy)).sqrt()).upper()
    norm = ((nominal.sqrt()+residual+reference_error)**2).upper()
    if not norm > 0:
        raise ArithmeticError('Missing nonzero physical family normalization.')
    payment = (arb(delta)*norm).upper()
    upper = (-arb(margin)+payment).upper()
    threshold = (arb(margin)/norm).lower()
    values = dict(mean_response_norm_squared=mean_norm,
        source_response_cross=cross, negative_mean_response_cross=negative_cross,
        weighted_response_11_upper=aa, weighted_response_22_upper=dd,
        weighted_response_abs12_upper=cc, response_two_by_two_eigenvalue_upper=top,
        negative_response_support_squared_upper=cross_energy,
        reference_nominal_family_norm_squared_upper=nominal,
        conditional_response_norm_upper=residual,
        inverse_reference_error_norm_upper=reference_error,
        directional_exact_family_norm_squared_upper=norm,
        complete_original_gain_upper=payment,
        complete_directional_family_quadratic_upper=upper,
        sufficient_original_gain_coefficient_lower=threshold)
    finite(values.values())
    return values


def consume(out):
    started = time.perf_counter()
    packets, guards = {}, {}
    def read(name, digest):
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted directional input: '+name)
        guards[name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
        return packets[name]
    family = read(FAMILY, FAMILY_SHA)
    for name in (PROFILE, CORE, CERT, METRIC, BUDGET, TRACE):
        read(name, family['packet_sha256'][name])
    profile, core, cert, metric, budget, trace = (packets[n] for n in (PROFILE, CORE, CERT, METRIC, BUDGET, TRACE))
    actual, projection = (read(n, trace['packet_sha256'][n]) for n in (ACTUAL, PROJECTION))
    if (family['schema'] != 'ACTUAL_C27_COHERENT_FAMILY_NORM_V1'
            or any(p.get('grid_bits', p.get('inherited_grid_bits')) != BITS for p in packets.values())
            or len(core['remaining_modes']) != 18 or cert['positive_dimension'] != 17
            or not cert['exact_reference_identities_checked'] or not profile['conditioning_certified']
            or any(p['remaining_modes'] != core['remaining_modes'] for p in packets.values())
            or any((p['source_integer_sha256'], p['source_bits'], p['source_row'])
                   != (SOURCE_ID, 384, 0) for p in (family, profile, metric, budget, trace))
            or profile['original_metadata_grid_bits'] != 512
            or profile['exact_far_allowance'] != trace['exact_far_allowance']
            or metric['positive_modes'] != profile['positive_modes'] or metric['basis'] != profile['basis']
            or any(trace['packet_sha256'][n] != guards[n] for n in (PROFILE, CORE, CERT, METRIC, BUDGET, ACTUAL, PROJECTION))
            or cert['packet_sha256']['core'] != guards[CORE]
            or any(p['input_sha256'] != core['input_sha256'] for p in (cert, actual, projection))
            or profile['profile_order'] != family['profile_order']):
        raise ValueError('Different same source, family, core, response, reference or normalization.')
    for name, digest in profile['inherited_dependency_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Different accepted helper: '+name)
        guards[name] = digest
    rational = lambda p, k: fmpq(*map(int, p['exact_rationals'][k]))
    scalar = lambda p, k: interval_ball(p['scalar_bounds'][k], BITS)
    eta = rational(budget, 'neumann_eta_upper')
    if (eta != rational(trace, 'neumann_eta_upper')
            or eta != fmpq(int(cert['scalar_bounds']['neumann_eta_upper'][1]), 1 << BITS)):
        raise ValueError('Different original two-sided inverse-reference allowance.')
    t_plus, fee, delta = (rational(family, k) for k in
        ('t_plus', 'actual_charge_floor', 'true_baseline_gain_coefficient'))
    if (t_plus != rational(budget, 'comparison_t_plus')
            or delta != fmpq(33, 16)*fmpq(profile['exact_far_allowance'])):
        raise ValueError('Different original optimistic coefficient or complete allowance.')
    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], BITS))
    reference = arb_mat([[arb(fmpq(v)) for v in row] for row in cert['reference_inverse']])
    euclidean = symmetric(matrix_in(profile['joint_euclidean_gram'], BITS))
    response, trace_q = directional_geometry(sigma, reference, euclidean,
        symmetric(matrix_in(metric['joint_inverse_metric'], BITS)), scalar(budget, 'metric_schur'),
        matrix_in([projection['matrices']['source_on_lift'][0]], BITS),
        matrix_in(actual['actual_finite_row_on_lift'], BITS),
        matrix_in(trace['reference_regression_actions'], BITS),
        symmetric(matrix_in(profile['conditional_gram'], BITS)))
    # Reuse, rather than recontract, the previously saved negative-column norm.
    response[0, 0] = scalar(family, 'negative_coordinate_euclidean_norm_squared')
    values = directional_budget(response, matrix_in(profile['measured_pair'], BITS),
        response[0, 0].upper(), scalar(family, 'exact_family_old_negative_margin_lower').lower(),
        scalar(family, 'exact_family_positive_energy_upper').upper(), t_plus, fee,
        rational(profile, 'remaining_energy_upper'), eta,
        scalar(family, 'reference_positive_euclidean_trace').upper(), max(arb(0), trace_q.upper()), delta)
    exported = {k: envelope(v, BITS) for k, v in values.items()}
    upper = fmpq(int(exported['complete_directional_family_quadratic_upper'][1]), 1 << BITS)
    for name, digest in guards.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Directional input changed during application: '+name)
    result = dict(schema='ACTUAL_C27_DIRECTIONAL_COHERENT_FAMILY_V1',
        packet_sha256={n:h for n,h in guards.items() if n in packets},
        inherited_dependency_sha256={n:h for n,h in guards.items() if n not in packets},
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_integer_sha256=SOURCE_ID, source_bits=384, source_row=0,
        remaining_modes=core['remaining_modes'], profile_order=profile['profile_order'],
        response_order=['negative_column', 'r_reference_response', 'bJ_reference_response',
                        'first_measured_reference_response', 'second_measured_reference_response'],
        response_gram=matrix_out(response, BITS), conditional_response_trace=envelope(trace_q, BITS),
        grid_bits=BITS, working_bits=ctx.prec, original_metadata_grid_bits=512,
        scalar_bounds=exported, exact_rationals={k:[str(v.p),str(v.q)] for k,v in dict(
            t_plus=t_plus, actual_charge_floor=fee, true_baseline_gain_coefficient=delta,
            complete_directional_family_upper=upper).items()},
        verdict=('same_stronger_base_family_excluded_by_directional_norm' if upper < 0
                 else 'inconclusive_directional_family_norm_budget'),
        new_response_gram_dimension=5, new_conditional_euclidean_trace=1,
        new_unsaved_reference_vector_products=2, reused_regression_reference_actions=2,
        new_inverse_solve_source_far_query=0, fixed_vector_selections=0,
        old_scientific_producers_or_fixtures_replayed=0,
        source_direction_parameter_precision_scan=False, full_entry_box_asserted_PSD=False,
        actual_retained_sign_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started, conditions=[
            'All inherited actual/common-source/domain,exact20frame,PSD physical and Euclidean geometry,conditional source,reference,complete original allowance and endpoint comparison premises remain.',
            'The SAME exact p=C^-1*k satisfies p*Lambda^-1*p<=p*C*p=k*C^-1*k<=Eplus. Neither p nor the actual Pplus inverse is replaced by a nominal reference solution.',
            'Good reference responses are e1 and e2+G^-1*m*. The18core responses use Q=T*R*T*. Five response columns retain all measured-source cross terms on the SAME H=A*A.',
            'The nominal reference family uses the exact source-dependent p. Its support uses the largest2x2response eigenvalue upper and the complete negative-column cross support.',
            'The SAME conditional residual q obeys q*q<=Delta*S. Its physical response norm is at most sqrt(Delta*tr(Hvv*Q*S*Q)), multiplied by the SAME |p2|<=sqrt(Eplus/fee). No independent source extrema are asserted jointly attainable.',
            'Only inverse-reference error uses tauP, multiplied by eta/(1-eta). R<=(1+eta)*Pplus^-1 and the exact positive energy imply ePR<=(1+eta)*Eplus. All norm cross terms are paid by triangle.',
            'A strict complete upper excludes only the SAME stronger endpoint comparison, uniformly over the inherited source set and dominated endpoint parameters. Nonnegative upper is inconclusive; no parameter,precision,source,direction or formula-refinement retry follows.',
            'This applies existing ellipse support,2x2spectral,Schur conditioning,relative inverse and norm inequalities. No generic novelty,maximal tolerance,actual retained/cofinal positivity,Robin,RH or Lean proof is claimed.'])
    write_packet(out, result)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Use a fresh output for one directional-family application.')
    ctx.prec = 1024
    consume(args.out)


if __name__ == '__main__':
    main()
