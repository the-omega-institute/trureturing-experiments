"""Bound the same coherent minimizer using its signed first-row equation.

All physical/reference Grams, energies and allowances are accepted inputs.
Only two mixed scalar dots and three small response contractions are new.
Schur elimination and convex quadratic endpoint bounds are standard methods.
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
DIRECTIONAL = 'full-gamma-directional-family-norm.json.zlib'
DIRECTIONAL_SHA = '6e12d77e100465978d3a9eb777280de15c1942256cbedf62c337fbe120f1c7d3'
FAMILY = 'full-gamma-coherent-family-norm.json.zlib'
PROFILE = 'full-gamma-conditioned-profile.json.zlib'
CORE = 'joint-gamma-6765-rho32-v2-schur.json'
CERT = 'full-gamma-directional-metric.json.zlib'
METRIC = 'full-gamma-finite-good-block-metric.json.zlib'
BUDGET = 'full-gamma-far-energy-budget.json.zlib'
TRACE = 'full-gamma-conditioned-joint-trace.json.zlib'
CAPACITY = 'full-gamma-negative-channel-capacity.json.zlib'
ACTUAL = 'full-gamma-actual-channel-projection.json.zlib'
PROJECTION = 'full-gamma-projection-upper.json.zlib'
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'


def finite(values):
    if any(not arb(v).is_finite() for v in values):
        raise ArithmeticError('Nonfinite exact-stationarity data.')


def hull(lo, hi):
    lo, hi = arb(lo), arb(hi)
    if lo > hi:
        raise ArithmeticError('Empty SAME-source stationarity interval.')
    return lo.union(hi)


def mixed_regression(sigma, rho, saved_actions):
    """Two NEW scalar dots; no reference vector action is replayed."""
    n = sigma.nrows()
    if (n < 2 or sigma.ncols() != n or (rho.nrows(), rho.ncols()) != (1, n)
            or (saved_actions.nrows(), saved_actions.ncols()) != (n-1, 2)):
        raise ValueError('Different stationarity frame dimensions.')
    finite(list(sigma.entries())+list(rho.entries())+list(saved_actions.entries()))
    if not sigma[0, 0] < 0:
        raise ArithmeticError('Missing inherited negative pivot.')
    rt = [rho[0, j+1]-rho[0, 0]*sigma[0, j+1]/sigma[0, 0] for j in range(n-1)]
    return arb_mat([[sum((rt[i]*saved_actions[i, j] for i in range(n-1)), arb(0))
                     for j in range(2)]])


def stationarity_budget(gram, measured, response, mixed, k1, k2, rr, rb,
        mean_energy, residual_energy, complete_energy, projection_energy,
        margin, energy, t_plus, fee, eta, c22_upper, delta, remaining,
        response_trace, inverse_error):
    if any((x.nrows(), x.ncols()) != dims for x, dims in
           [(gram, (2, 2)), (measured, (1, 2)), (response, (5, 5)), (mixed, (1, 2))]):
        raise ValueError('Different small exact-stationarity dimensions.')
    scalars = (k1, k2, rr, rb, mean_energy, residual_energy, complete_energy,
        projection_energy, margin, energy, t_plus, fee, eta, c22_upper,
        delta, remaining, response_trace, inverse_error)
    finite([v for x in (gram, measured, response, mixed) for v in x.entries()]+list(scalars))
    if (not gram[0, 0] > 0 or not margin > 0 or not t_plus > 0 or not fee > 0
            or not 0 <= eta < 1 or not mean_energy >= 0
            or any(not x >= 0 for x in (rr, residual_energy, complete_energy,
                  projection_energy, energy, delta, remaining, response_trace, inverse_error))):
        raise ArithmeticError('Missing SAME-family PSD, energy or allowance domain.')
    inverse_rr_error = (arb(eta)/(1-arb(eta))*rr.upper()).upper()
    aa = 1/arb(t_plus)+gram[0, 0]+rr+arb(0, inverse_rr_error)
    cross_center = gram[0, 1]+measured[0, 0]+rb+(mixed*measured.transpose())[0, 0]
    conditional_cross = (rr.upper()*arb(residual_energy)).sqrt().upper()
    inverse_cross = (arb(eta)/(1-arb(eta))
                     *(rr.upper()*arb(complete_energy)).sqrt()).upper()
    cross = cross_center+arb(0, conditional_cross+inverse_cross)
    good_second = arb(fee)+gram[1, 1]+2*measured[0, 1]+projection_energy
    mean_min = max(arb(0), mean_energy.lower().sqrt()-arb(residual_energy).sqrt()).lower()
    second_lower = (good_second.lower()+mean_min**2/(1+arb(eta))).lower()
    second = hull(second_lower, c22_upper)
    ratio, base = cross/aa, k1/aa
    schur = second-cross**2/aa
    schur = hull(max(arb(fee), schur.lower()), schur.upper())
    if not aa > 0 or not schur > 0:
        raise ArithmeticError('No certified positive stationarity denominator.')
    numerator = k2-ratio*k1
    p2_ratio = numerator/schur
    energy_remainder = (arb(energy)-(k1**2/aa).lower()).upper()
    if not energy_remainder >= 0:
        raise ArithmeticError('Exact first-row energy contradicts inherited family energy.')
    p2_radius = (energy_remainder/arb(fee)).sqrt().upper()
    p2 = hull(max(p2_ratio.lower(), -p2_radius), min(p2_ratio.upper(), p2_radius))
    p1 = base-ratio*p2
    h0 = arb_mat([[1], [-base], [0], [0], [0]])
    h2 = arb_mat([[0], [ratio], [-1], [-measured[0, 0]], [-measured[0, 1]]])
    # Three NEW contractions of the saved Gram. The exact norm is convex in p2.
    a0 = (h0.transpose()*response*h0)[0, 0]
    b0 = (h0.transpose()*response*h2)[0, 0]
    c0 = (h2.transpose()*response*h2)[0, 0]
    endpoints = [(a0+2*b0*p+c0*p**2).upper() for p in (p2.lower(), p2.upper())]
    nominal = max(endpoints)
    if nominal < 0 or c0.upper() < 0:
        raise ArithmeticError('Saved exact PSD response Gram has inconsistent upper.')
    residual = (abs(p2).upper()
                *(arb(remaining)*arb(response_trace)).sqrt()).upper()
    norm = ((nominal.sqrt()+residual+arb(inverse_error))**2).upper()
    if not norm > 0:
        raise ArithmeticError('Missing nonzero exact-family normalization.')
    payment = (arb(delta)*norm).upper()
    values = dict(first_source_diagonal=aa, source_cross_center=cross_center,
        conditional_source_cross_radius_upper=conditional_cross,
        inverse_source_cross_radius_upper=inverse_cross, complete_source_cross=cross,
        second_source_diagonal=second, source_cross_ratio=ratio,
        first_source_base_coefficient=base, stationarity_schur=schur,
        stationarity_second_numerator=numerator, second_coefficient_ratio_interval=p2_ratio,
        first_row_energy_remainder_upper=energy_remainder,
        energy_second_coefficient_radius_upper=p2_radius,
        exact_second_coefficient=p2, exact_first_coefficient=p1,
        nominal_quadratic_constant=a0, nominal_quadratic_half_linear=b0,
        nominal_quadratic_square=c0, nominal_lower_endpoint_upper=endpoints[0],
        nominal_upper_endpoint_upper=endpoints[1],
        stationarity_nominal_norm_squared_upper=nominal,
        stationarity_conditional_response_norm_upper=residual,
        inherited_inverse_reference_error_norm_upper=arb(inverse_error),
        complete_stationarity_norm_squared_upper=norm,
        complete_original_gain_upper=payment,
        complete_stationarity_family_quadratic_upper=(-arb(margin)+payment).upper(),
        sufficient_original_gain_coefficient_lower=(arb(margin)/norm).lower())
    finite(values.values())
    return values


def consume(out):
    started = time.perf_counter()
    packets, guards = {}, {}
    def read(name, digest):
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted stationarity input: '+name)
        guards[name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
        return packets[name]
    directional = read(DIRECTIONAL, DIRECTIONAL_SHA)
    family = read(FAMILY, directional['packet_sha256'][FAMILY])
    for name in (PROFILE, CORE, CERT, METRIC, BUDGET, TRACE, ACTUAL, PROJECTION):
        read(name, directional['packet_sha256'][name])
    capacity = read(CAPACITY, family['packet_sha256'][CAPACITY])
    profile, core, cert, metric, budget, trace, actual, projection = (
        packets[n] for n in (PROFILE, CORE, CERT, METRIC, BUDGET, TRACE, ACTUAL, PROJECTION))
    if (directional['schema'] != 'ACTUAL_C27_DIRECTIONAL_COHERENT_FAMILY_V1'
            or family['schema'] != 'ACTUAL_C27_COHERENT_FAMILY_NORM_V1'
            or directional['response_order'] != ['negative_column', 'r_reference_response',
                'bJ_reference_response', 'first_measured_reference_response', 'second_measured_reference_response']
            or any(p.get('grid_bits', p.get('inherited_grid_bits')) != BITS
                   or p['remaining_modes'] != core['remaining_modes'] for p in packets.values())
            or len(core['remaining_modes']) != 18 or cert['positive_dimension'] != 17
            or not cert['exact_reference_identities_checked'] or not profile['conditioning_certified']
            or profile['original_metadata_grid_bits'] != 512
            or any((p['source_integer_sha256'], p['source_bits'], p['source_row'])
                   != (SOURCE_ID, 384, 0) for p in
                   (directional, family, profile, metric, budget, trace, capacity))
            or any(p['profile_order'] != profile['profile_order'] for p in (directional, family))
            or profile['exact_far_allowance'] != trace['exact_far_allowance']
            or metric['basis'] != profile['basis'] or metric['positive_modes'] != profile['positive_modes']
            or any(trace['packet_sha256'][n] != guards[n] for n in
                   (PROFILE, CORE, CERT, METRIC, BUDGET, ACTUAL, PROJECTION, CAPACITY))
            or cert['packet_sha256']['core'] != guards[CORE]
            or any(p['input_sha256'] != core['input_sha256'] for p in (cert, actual, projection))):
        raise ValueError('Different SAME source, frame, conditional geometry or original grid.')
    for name, digest in profile['inherited_dependency_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Different inherited helper: '+name)
        guards[name] = digest
    rational = lambda p, k: fmpq(*map(int, p['exact_rationals'][k]))
    scalar = lambda p, k: interval_ball(p['scalar_bounds'][k], BITS)
    eta = rational(budget, 'neumann_eta_upper')
    t_plus, fee, delta = (rational(family, k) for k in
        ('t_plus', 'actual_charge_floor', 'true_baseline_gain_coefficient'))
    if (eta != rational(trace, 'neumann_eta_upper')
            or eta != fmpq(int(cert['scalar_bounds']['neumann_eta_upper'][1]), 1 << BITS)
            or t_plus != rational(budget, 'comparison_t_plus')
            or delta != fmpq(33, 16)*fmpq(profile['exact_far_allowance'])):
        raise ValueError('Different original endpoint, reference or complete allowance.')
    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], BITS))
    rho = matrix_in([projection['matrices']['source_on_lift'][0]], BITS)
    mixed = mixed_regression(sigma, rho, matrix_in(trace['reference_regression_actions'], BITS))
    k2 = (matrix_in(actual['actual_finite_row_on_lift'], BITS)[0, 0]
          +matrix_in(profile['intersected_far_row'], BITS)[0, 0])
    values = stationarity_budget(symmetric(matrix_in(metric['joint_inverse_metric'], BITS)),
        matrix_in(profile['measured_pair'], BITS), matrix_in(directional['response_gram'], BITS),
        mixed, rho[0, 0], k2, scalar(capacity, 'steering_reference_rRr'),
        scalar(capacity, 'steering_reference_rRb'), scalar(trace, 'mean_reference_energy'),
        scalar(trace, 'conditional_residual_reference_energy_upper').upper(),
        scalar(trace, 'complete_reference_energy_upper').upper(),
        scalar(budget, 'far_energy_projection_floor'),
        scalar(family, 'exact_family_old_negative_margin_lower').lower(),
        scalar(family, 'exact_family_positive_energy_upper').upper(), t_plus, fee, eta,
        scalar(family, 'positive_source_middle_22_upper').upper(), delta,
        rational(profile, 'remaining_energy_upper'),
        max(arb(0), interval_ball(directional['conditional_response_trace'], BITS).upper()),
        scalar(directional, 'inverse_reference_error_norm_upper').upper())
    exported = {k: envelope(v, BITS) for k, v in values.items()}
    upper = fmpq(int(exported['complete_stationarity_family_quadratic_upper'][1]), 1 << BITS)
    for name, digest in guards.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Accepted input changed during stationarity application: '+name)
    result = dict(schema='ACTUAL_C27_STATIONARITY_COHERENT_FAMILY_V1',
        packet_sha256={n:h for n,h in guards.items() if n in packets},
        inherited_dependency_sha256={n:h for n,h in guards.items() if n not in packets},
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_integer_sha256=SOURCE_ID, source_bits=384, source_row=0,
        remaining_modes=core['remaining_modes'], profile_order=profile['profile_order'],
        response_order=directional['response_order'], mixed_regression=matrix_out(mixed, BITS),
        grid_bits=BITS, working_bits=ctx.prec, original_metadata_grid_bits=512,
        scalar_bounds=exported, exact_rationals={k:[str(v.p),str(v.q)] for k,v in dict(
            t_plus=t_plus, actual_charge_floor=fee, true_baseline_gain_coefficient=delta,
            complete_stationarity_family_upper=upper).items()},
        verdict=('same_stronger_base_family_excluded_by_stationarity'
                 if upper < 0 else 'inconclusive_stationarity_family_norm_budget'),
        new_mixed_scalar_dots=2, new_saved_gram_scalar_contractions=3,
        new_Gram_trace_reference_vector_action_inverse_solve_source_far_query=0,
        fixed_vector_selections=0, old_scientific_producers_or_fixtures_replayed=0,
        source_direction_parameter_precision_scan=False, full_entry_box_asserted_PSD=False,
        actual_retained_sign_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started, conditions=[
            'ALL inherited actual/common-source,domain,exact20frame,physical and Euclidean PSD,conditional source,two-sided relative inverse,complete original512metadata allowance and endpoint comparison premises remain.',
            'The SAME exact p=C^-1*k and v=-C0^-1*F**p are retained. Signed first-row stationarity and the positive Schur complement enclose p2 and p1,not a nominal fixed vector or a free coefficient ellipse.',
            'c,k2,mean response,m,q belong to the SAME source. Interval hulls safely enclose their dependence; no independently extremal interval corners are asserted jointly attainable.',
            'Reference energies, reverse triangle and the two-sided factor give the second diagonal lower; C11 uses the saved rr and eta. Cauchy bounds only the conditional cross and inverse-reference error.',
            'For each SAME realization the saved response-Gram quadratic is convex in REAL p2. Its two endpoint uppers bound the nominal norm even though x,c,m and p2 are dependent. Full entry boxes are not asserted PSD.',
            'The saved conditional Euclidean trace pays maxabs(p2)*sqrt(Delta*tauQ). The prior inverse-reference error bound still applies to this exact family energy. All norm cross terms are paid by triangle.',
            'Strict complete upper excludes only the SAME stronger endpoint comparison,uniformly over inherited source set and dominated endpoint parameters. A nonnegative upper is inconclusive. One application,no parameter/source/direction/precision or formula retry.',
            'Standard Schur elimination,PSD quadratic convexity,energy and norm inequalities are reused. This source-specific application claims no generic novelty,maximal tolerance,actual retained/cofinal positivity,Robin,RH or Lean proof.'])
    write_packet(out, result)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Use one fresh .zlib output; no stationarity-budget retry.')
    ctx.prec = 1024
    consume(args.out)


if __name__ == '__main__':
    main()
