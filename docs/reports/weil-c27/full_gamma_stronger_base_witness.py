"""Certify one fixed twenty-profile witness after paying the true-base gain.

Reuse every old source, solve and Gram. A new combined reference action
selects exact dyadic coefficients; the original enclosing forms, not the
selection midpoint or reference inverse, certify their quadratic value.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import zlib

from flint import arb, arb_mat, ctx, fmpq
from diagonal_floor_witness import interval_ball
from joint_gamma_boundary import matrix_in, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
PINS = {
    'full-gamma-conditioned-profile.json.zlib': '2e1c11ef3807445ce6b157d47232bc6086ac96f8a07786e4b03b7355e61c00ae',
    'joint-gamma-6765-rho32-v2-schur.json': '13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd',
    'full-gamma-directional-metric.json.zlib': 'f6bea1b9b43c5ca8c85fc6cd7222038f2a40b81ed1ce2ef873af91ce1ca987a2',
    'full-gamma-actual-channel-projection.json.zlib': '63e3c2e291940ca3549fa4a1303b55b4d92c08d594a50daf7ca2132087ad2745',
    'full-gamma-projection-upper.json.zlib': 'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3',
    'full-gamma-finite-good-block-metric.json.zlib': '9ba53b70a0f1cad25b9e92894eff93f1340dc31d1a9bf945c560b6dc6940d072',
    'full-gamma-far-energy-budget.json.zlib': '11ccabd6110436b2b29f1a910c6c36a5fd8479af19d183cf614826ed9ad746c6',
    'full-gamma-actual-joint-lower.json.zlib': 'f294428dc9685e04fca777a5dc71ac0374fc17a4c69246422e8b4325351da5a9',
    'full-gamma-conditioned-joint-trace.json.zlib': '6656a9e9ab05dd090168c93a359eeb8930998b89e855755b062f2493634dd0dd',
    'full-gamma-negative-channel-capacity.json.zlib': '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
}
PROFILE, CORE, CERT, ACTUAL, PROJECTION, METRIC, BUDGET, JOINT, TRACE, CAPACITY = PINS


def point(value):
    lo, hi = envelope(value, BITS)
    return fmpq(int(lo)+int(hi), 1 << (BITS+1))


def points(matrix):
    return arb_mat([[arb(point(matrix[i, j])) for j in range(matrix.ncols())]
                    for i in range(matrix.nrows())])


def inverse2(matrix):
    """Selection only: explicit inverse of a positive definite2x2 point matrix."""
    det = matrix[0, 0]*matrix[1, 1]-matrix[0, 1]*matrix[1, 0]
    if (matrix.nrows(), matrix.ncols()) != (2, 2) or not matrix[0, 0] > 0 or not det > 0:
        raise ArithmeticError('Missing positive two-dimensional selection denominator.')
    return arb_mat([[matrix[1, 1], -matrix[0, 1]],
                    [-matrix[1, 0], matrix[0, 0]]])/det


def witness_upper(w, gram, sigma, euclidean, rho, beta, regression,
                  measured, conditional, remaining, t_plus, fee, delta):
    n = sigma.nrows()
    if ((w.nrows(), w.ncols()) != (n+2, 1)
            or (sigma.nrows(), sigma.ncols()) != (n, n)
            or (gram.nrows(), gram.ncols()) != (2, 2)
            or (euclidean.nrows(), euclidean.ncols()) != (n+2, n+2)
            or (rho.nrows(), rho.ncols()) != (1, n)
            or (beta.nrows(), beta.ncols()) != (1, n)
            or (regression.nrows(), regression.ncols()) != (2, n)
            or (measured.nrows(), measured.ncols()) != (1, 2)
            or (conditional.nrows(), conditional.ncols()) != (n, n)
            or not remaining >= 0 or not t_plus > 0 or not fee > 0 or not delta >= 0):
        raise ValueError('Incomplete fixed-witness source, geometry or scalar domain.')
    g = arb_mat([[w[i, 0]] for i in range(2)])
    v = arb_mat([[w[i+2, 0]] for i in range(n)])
    base = (g.transpose()*gram*g+v.transpose()*sigma*v)[0, 0]
    norm_sq = (w.transpose()*euclidean*w)[0, 0]
    variance = (v.transpose()*conditional*v)[0, 0]
    if norm_sq.upper() < 0 or variance.upper() < 0:
        raise ArithmeticError('The inherited exact PSD forms have inconsistent enclosures.')
    r_action = gram[0, 0]*g[0, 0]+gram[0, 1]*g[1, 0]+(rho*v)[0, 0]
    # Collect the SAME measured pair before enclosing its contribution.
    b_mean = (gram[1, 0]*g[0, 0]+gram[1, 1]*g[1, 0]
              +(beta*v)[0, 0]+(measured*(g+regression*v))[0, 0])
    residual_cap = (arb(remaining)*max(arb(0), variance.upper())).sqrt().upper()
    b_cap = (abs(b_mean).upper()+residual_cap).upper()
    gap = (arb(delta)*max(arb(0), norm_sq.upper())).upper()
    r_fee = (arb(t_plus)*abs(r_action).upper()**2).upper()
    b_fee = (b_cap**2/arb(fee)).upper()
    upper = (base.upper()+gap+r_fee+b_fee).upper()
    result = dict(old_base_quadratic=base, euclidean_norm_squared=norm_sq,
        conditional_variance=variance, complete_r_action=r_action,
        conditional_b_mean=b_mean, conditional_b_residual_cap=residual_cap,
        complete_b_action_magnitude_upper=b_cap, true_baseline_gain_upper=gap,
        shift_update_upper=r_fee, charge_update_upper=b_fee,
        complete_stronger_target_quadratic_upper=upper)
    if any(not x.is_finite() for x in result.values()):
        raise ArithmeticError('Nonfinite directed fixed-witness result.')
    return result


def consume(out):
    started = time.perf_counter()
    packets, guarded = {}, {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Saved input hash mismatch: '+name)
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
        guarded[name] = digest
    profile, core, cert, actual, projection, metric, budget, joint, trace, capacity = packets.values()
    if (any(p['grid_bits'] != BITS or p['remaining_modes'] != core['remaining_modes']
            for p in packets.values()) or len(core['remaining_modes']) != 18
            or not profile['conditioning_certified'] or not cert['exact_reference_identities_checked']
            or profile['original_metadata_grid_bits'] != 512
            or any((p['source_integer_sha256'], p['source_bits'], p['source_row'])
                   != (SOURCE_ID, 384, 0) for p in (profile, metric, budget, trace, capacity))):
        raise ValueError('Different source, core, exact-reference or conditioning roles.')
    for name in (PROFILE, CORE, CERT, ACTUAL, PROJECTION, METRIC, BUDGET, CAPACITY):
        if trace['packet_sha256'][name] != PINS[name]:
            raise ValueError('Different accepted joint trace dependency: '+name)
    if (profile['packet_sha256'][METRIC] != PINS[METRIC]
            or joint['packet_sha256']['actual'] != PINS[ACTUAL]
            or joint['packet_sha256']['core'] != PINS[CORE]
            or profile['exact_far_allowance'] != trace['exact_far_allowance']
            or metric['positive_modes'] != profile['positive_modes']
            or metric['basis'] != profile['basis']):
        raise ValueError('Different common charge, full allowance, frame or Gram.')
    for name, digest in profile['inherited_dependency_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Inherited helper hash mismatch: '+name)
        guarded[name] = digest
    schemas = ['ACTUAL_C27_CONDITIONED_TWENTY_PROFILE_GRAM_V1',
        'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1', 'ACTUAL_C27_FULL_GAMMA_DIRECTIONAL_METRIC_V1',
        'ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1', 'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1',
        'ACTUAL_C27_FINITE_GOOD_BLOCK_JOINT_METRIC_V1', 'ACTUAL_C27_FAR_INVERSE_ENERGY_BUDGET_V1',
        'ACTUAL_C27_FULL_GAMMA_ACTUAL_JOINT_LOWER_V1', 'ACTUAL_C27_CONDITIONED_JOINT_TRACE_BUDGET_V1',
        'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1']
    if [p['schema'] for p in packets.values()] != schemas:
        raise ValueError('Different accepted packet schemas.')
    rational = lambda p, k: fmpq(*map(int, p['exact_rationals'][k]))
    scalar = lambda p, k: interval_ball(p['scalar_bounds'][k], BITS)
    t_plus = rational(budget, 'comparison_t_plus')
    fee = fmpq(int(joint['scalar_bounds']['actual_common_charge_floor'][0]), 1 << BITS)
    remaining = rational(profile, 'remaining_energy_upper')
    delta = fmpq(33, 16)*fmpq(profile['exact_far_allowance'])
    sigma = matrix_in(core['reduced_schur_intervals'], BITS)
    gram = matrix_in(metric['joint_inverse_metric'], BITS)
    rho = matrix_in([projection['matrices']['source_on_lift'][0]], BITS)
    beta = matrix_in(actual['actual_finite_row_on_lift'], BITS)
    regression, measured = matrix_in(profile['regression'], BITS), matrix_in(profile['measured_pair'], BITS)
    conditional = matrix_in(profile['conditional_gram'], BITS)
    euclidean = matrix_in(profile['joint_euclidean_gram'], BITS)
    pivot = sigma[0, 0]
    if not pivot < 0 or not fee > 0 or not delta > 0:
        raise ArithmeticError('Missing inherited pivot, actual charge floor or true-base allowance.')
    transport = arb_mat([[-sigma[0, j+1]/pivot if i == 0 else arb(int(i == j+1))
                          for j in range(17)] for i in range(18)])
    reference = arb_mat([[arb(fmpq(v)) for v in row] for row in cert['reference_inverse']])
    # Select one nominal Woodbury correction using saved positive energies.
    cross = scalar(capacity, 'steering_reference_rRb')
    cross += ((rho*transport)*matrix_in(trace['reference_regression_actions'], BITS)
              *measured.transpose())[0, 0]
    inverse_good = inverse2(points(gram))
    observed = points(measured)
    good_rows = arb_mat([[gram[0, 0], gram[0, 1]],
                         [gram[1, 0]+measured[0, 0], gram[1, 1]+measured[0, 1]]])
    q_projection = (observed*inverse_good*observed.transpose())[0, 0]
    middle = points(arb_mat([[1/arb(t_plus)+gram[0, 0]+scalar(capacity, 'steering_reference_rRr'),
                             gram[0, 1]+measured[0, 0]+cross],
                            [gram[0, 1]+measured[0, 0]+cross,
                             arb(fee)+gram[1, 1]+2*measured[0, 1]+q_projection
                             +scalar(trace, 'mean_reference_energy')]]))
    mean = beta+measured*regression
    negative_row = points(arb_mat([[rho[0, 0]], [mean[0, 0]]]))
    coefficient = inverse2(middle)*negative_row
    good_coefficients = -inverse_good*points(good_rows).transpose()*coefficient
    new_rhs = points((rho*transport).transpose())*coefficient[0, 0]
    new_rhs += points((mean*transport).transpose())*coefficient[1, 0]
    positive_coefficients = -reference*new_rhs
    lifted_coefficients = points(transport)*positive_coefficients
    lifted_coefficients[0, 0] += 1
    raw_coefficients = [good_coefficients[i, 0] for i in range(2)]
    raw_coefficients += [lifted_coefficients[i, 0] for i in range(18)]
    exact_w = [point(v) for v in raw_coefficients]
    w = arb_mat([[arb(v)] for v in exact_w])
    values = witness_upper(w, gram, sigma, euclidean, rho, beta, regression,
        measured, conditional, remaining, t_plus, fee, delta)
    pairs = {k: envelope(v, BITS) for k, v in values.items()}
    upper = fmpq(int(pairs['complete_stronger_target_quadratic_upper'][1]), 1 << BITS)
    for name, digest in guarded.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Material input changed during the new witness application: '+name)
    result = dict(schema='ACTUAL_C27_STRONGER_BASE_FIXED_WITNESS_V1', packet_sha256=PINS,
        inherited_dependency_sha256=profile['inherited_dependency_sha256'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_integer_sha256=SOURCE_ID, source_bits=384, source_row=0,
        remaining_modes=core['remaining_modes'], profile_order=profile['profile_order'],
        grid_bits=BITS, working_bits=ctx.prec, original_metadata_grid_bits=512,
        fixed_exact_dyadic_coefficients=[str(v) for v in exact_w], scalar_bounds=pairs,
        exact_rationals={k: [str(v.p), str(v.q)] for k, v in dict(t_plus=t_plus,
            actual_charge_floor=fee, remaining_energy_upper=remaining,
            true_baseline_gain_coefficient=delta, complete_quadratic_upper=upper).items()},
        verdict=('same_stronger_base_source_family_excluded_by_fixed_witness'
                 if upper < 0 else 'inconclusive_fixed_stronger_base_witness'),
        new_combined_reference_actions=1, new_selection_two_by_two_inverses=2,
        new_source_far_target_inverse_solve_queries=0,
        old_scientific_producers_or_fixtures_replayed=0,
        source_direction_parameter_precision_cutoff_scan=False,
        full_entry_box_asserted_PSD=False, midpoint_reference_used_for_certification=False,
        actual_retained_form_sign_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started,
        conditions=[
            'All inherited actual common-form/domain, SAME source/lift, normalization, physical-majorant, conditional-residual, source-energy and complete-far allowance premises remain.',
            'A=[iota_g Y,V] is the SAME saved20profile frame, with A*Lstar*A=diag(GJJ,Sigma), and its saved Euclidean Gram is A*A. The full entry box is not asserted PSD.',
            'The SAME true stronger baseline obeys 0<=Lbar-Lstar<=delta I, delta=2*kappa*epsilon with the complete original allowance; the two padding grids are not identified.',
            'The true rows on A are rA=[a,c,rho] and bA=[c+u,d+v,betaJ+y], with measured m=(u,v) and y=mB+q, qq*<=remaining*S. The observed contributions are collected before enclosure.',
            'Selection uses nominal saved energies and one new combined reference action only to choose a fixed exact dyadic coefficient vector. The original interval base and Euclidean forms and SAME conditional ellipsoid certify that vector; no inverse error payment or midpoint equality is assumed.',
            'A strict negative full quadratic upper gives a nonzero negative direction for the SAME stronger lower construction at optimistic t_plus,Dlo. Loewner order extends this to0<t<=t_plus,D>=K>=Dlo, including the supplied eligible shifts; endpoint eligibility is unnecessary.',
            'A nonnegative witness upper is inconclusive and triggers no direction/parameter/precision retry. Other sources, other majorants, larger shifts, actual retained/cofinal positivity or negativity, Robin and RH remain unresolved. No generic novelty or Lean-proof claim.'])
    write_packet(out, result)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Use a fresh .zlib output; no direction or precision retry.')
    ctx.prec = 1024
    consume(args.out)


if __name__ == '__main__':
    main()
