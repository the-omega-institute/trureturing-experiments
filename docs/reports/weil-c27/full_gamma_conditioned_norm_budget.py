"""Consume saved reference energies to bound the complete joint comparison.

PSD trace monotonicity and triangle inequalities replace the unresolved
joint-gradient contraction. All paid center and inverse-error bounds are
reused; no vector/matrix action, inverse, solve or source query runs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import zlib

from flint import arb, ctx, fmpq
from diagonal_floor_witness import interval_ball
from joint_gamma_boundary import write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
PINS = {
    'full-gamma-conditioned-joint-trace.json.zlib': '6656a9e9ab05dd090168c93a359eeb8930998b89e855755b062f2493634dd0dd',
    'full-gamma-conditioned-profile.json.zlib': '2e1c11ef3807445ce6b157d47232bc6086ac96f8a07786e4b03b7355e61c00ae',
    'joint-gamma-6765-rho32-v2-schur.json': '13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd',
    'full-gamma-negative-channel-capacity.json.zlib': '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
    'full-gamma-actual-channel-projection.json.zlib': '63e3c2e291940ca3549fa4a1303b55b4d92c08d594a50daf7ca2132087ad2745',
    'full-gamma-projection-upper.json.zlib': 'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3',
    'full-gamma-finite-good-block-metric.json.zlib': '9ba53b70a0f1cad25b9e92894eff93f1340dc31d1a9bf945c560b6dc6940d072',
    'full-gamma-far-energy-budget.json.zlib': '11ccabd6110436b2b29f1a910c6c36a5fd8479af19d183cf614826ed9ad746c6',
}
TRACE, PROFILE, CORE, CAPACITY, ACTUAL, PROJECTION, METRIC, BUDGET = PINS


def norm_budget(pivot, s00, rho0, mean0, observed_offdiagonal,
                mean_reference, rho_reference, conditional_trace, remaining,
                nominal_a, center, complete_fee):
    """A uniform complete-budget upper from scalar PSD reference bounds.

    Write Ghat=W-d*e0*e0^T, d=-1/pivot>0. Under W,S PSD,
    W*S*W <= tr(S*W)*W, even when either matrix is singular.
    Hence ||S^(1/2)Ghat*v|| <= sqrt(tr(SW)*vWv)+d*|v0|sqrt(S00).
    Dropping -rWr in the nominal curvature is an upper bound. The saved
    complete_fee already prices actual-minus-nominal on the same source.
    """
    inputs = [pivot, s00, rho0, mean0, observed_offdiagonal,
              mean_reference, rho_reference, conditional_trace, arb(remaining),
              arb(nominal_a), center, complete_fee]
    if any(not x.is_finite() for x in inputs):
        raise ArithmeticError('Nonfinite saved scalar norm-budget input.')
    if (not pivot < 0 or not nominal_a > 0 or not remaining >= 0
            or any(x.upper() < 0 for x in
                   (s00, mean_reference, rho_reference, conditional_trace, complete_fee))):
        raise ArithmeticError('Missing inherited PSD, pivot, denominator or budget domain.')
    d = (1/(-pivot)).upper()
    s = max(arb(0), s00.upper())
    mean = max(arb(0), mean_reference.upper())
    rho = max(arb(0), rho_reference.upper())
    trace = max(arb(0), conditional_trace.upper())
    rho_first, mean_first = abs(rho0).upper(), abs(mean0).upper()
    rho_norm = ((trace*rho).sqrt()+d*rho_first*s.sqrt()).upper()
    mean_norm = ((trace*mean).sqrt()+d*mean_first*s.sqrt()).upper()
    hcap = (abs(observed_offdiagonal).upper()
            +d*rho_first*mean_first+(rho*mean).sqrt()).upper()
    linear = (arb(remaining).sqrt()*(hcap*rho_norm/arb(nominal_a)+mean_norm)).upper()
    curvature = (rho_norm**2/arb(nominal_a)+d*s).upper()
    upper = center.upper()+2*linear+arb(remaining)*curvature+complete_fee.upper()
    values = dict(negative_pivot_magnitude_inverse_upper=d,
        conditional_first_coordinate_variance_upper=s,
        mean_first_coordinate=mean0, source_first_coordinate=rho0,
        observed_offdiagonal=observed_offdiagonal,
        reused_mean_reference_energy=mean_reference,
        reused_source_reference_energy=rho_reference,
        reused_conditional_reference_trace=conditional_trace,
        source_conditional_norm_upper=rho_norm, mean_conditional_norm_upper=mean_norm,
        nominal_offdiagonal_magnitude_upper=hcap, linear_support_upper=linear,
        curvature_upper=curvature, reused_nominal_center=center,
        reused_complete_reference_perturbation_upper=complete_fee,
        complete_joint_budget_upper=upper)
    if any(not x.is_finite() for x in values.values()):
        raise ArithmeticError('Nonfinite complete scalar norm-budget output.')
    return values


def consume(out):
    started = time.perf_counter()
    packets = {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Input hash mismatch: '+name)
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    trace, profile, core, capacity, actual, projection, metric, budget = packets.values()
    if [p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_CONDITIONED_JOINT_TRACE_BUDGET_V1',
            'ACTUAL_C27_CONDITIONED_TWENTY_PROFILE_GRAM_V1',
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1',
            'ACTUAL_C27_FIXED_NEGATIVE_CHANNEL_CAPACITY_V1',
            'ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
            'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1',
            'ACTUAL_C27_FINITE_GOOD_BLOCK_JOINT_METRIC_V1',
            'ACTUAL_C27_FAR_INVERSE_ENERGY_BUDGET_V1']:
        raise ValueError('Different scalar comparison packet roles.')
    if (any(trace['packet_sha256'][name] != PINS[name] for name in list(PINS)[1:])
            or any(p['grid_bits'] != BITS for p in packets.values())
            or any(p['remaining_modes'] != core['remaining_modes']
                   for p in (trace, profile, capacity, actual, projection, metric, budget))
            or len(core['remaining_modes']) != 18
            or trace['original_metadata_grid_bits'] != 512
            or any((p['source_integer_sha256'], p['source_bits'], p['source_row'])
                   != (SOURCE_ID, 384, 0) for p in (trace, profile, capacity, metric, budget))
            or not profile['conditioning_certified']):
        raise ValueError('Different inherited comparison, grids, conditioning or source.')
    helpers = profile['inherited_dependency_sha256']
    for name, digest in helpers.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Inherited helper hash mismatch: '+name)
    def scalar(pair):
        return interval_ball(pair, BITS)
    scalars = trace['scalar_bounds']
    measured = [scalar(v) for v in profile['measured_pair'][0]]
    mean0 = (scalar(actual['actual_finite_row_on_lift'][0][0])
             +sum((measured[j]*scalar(profile['regression'][j][0]) for j in range(2)), arb(0)))
    remaining = fmpq(*map(int, trace['exact_rationals']['remaining_energy_upper']))
    nominal_a = fmpq(*map(int, trace['exact_rationals']['nominal_border_denominator']))
    values = norm_budget(scalar(core['reduced_schur_intervals'][0][0]),
        scalar(profile['conditional_gram'][0][0]),
        scalar(projection['matrices']['source_on_lift'][0][0]), mean0,
        scalar(metric['joint_inverse_metric'][0][1])+measured[0],
        scalar(scalars['mean_reference_energy']),
        scalar(capacity['scalar_bounds']['steering_reference_rRr']),
        scalar(scalars['conditional_reference_trace']), remaining, nominal_a,
        scalar(scalars['nominal_center_quadratic']),
        scalar(scalars['complete_reference_perturbation_upper']))
    pairs = {k: envelope(v, BITS) for k, v in values.items()}
    upper = fmpq(int(pairs['complete_joint_budget_upper'][1]), 1 << BITS)
    old_upper = fmpq(*map(int, trace['exact_rationals']['complete_joint_budget_upper']))
    a_pair = budget['scalar_bounds']['border_a']
    a_lower = fmpq(int(a_pair[0]), 1 << BITS)
    if not a_lower > 0:
        raise ArithmeticError('The inherited actual comparison denominator must be positive.')
    exact = dict(remaining_energy_upper=remaining, nominal_border_denominator=nominal_a,
        original_joint_trace_upper=old_upper, complete_joint_budget_upper=upper,
        actual_comparison_denominator_lower=a_lower)
    if upper < 0:
        exact['fixed_comparison_determinant_positive_lower'] = a_lower*(-upper)
    result = dict(schema='ACTUAL_C27_CONDITIONED_SCALAR_NORM_BUDGET_V1',
        packet_sha256=PINS, inherited_dependency_sha256=helpers,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_integer_sha256=SOURCE_ID, source_bits=384, source_row=0,
        remaining_modes=core['remaining_modes'], grid_bits=BITS, working_bits=ctx.prec,
        original_metadata_grid_bits=512,
        scalar_bounds=pairs, exact_rationals={k: [str(v.p), str(v.q)] for k, v in exact.items()},
        verdict=('fixed_optimistic_comparison_excluded_by_scalar_norm'
                 if upper < 0 else 'inconclusive_conditional_scalar_norm_budget'),
        strictly_tighter_complete_budget_upper=bool(upper < old_upper),
        comparison_coefficient_eligibility_claim=False,
        new_scalar_first_coordinate_mean=True,
        new_source_far_target_inverse_solve_queries=0,
        new_vector_matrix_reference_actions_or_contractions=0,
        old_scientific_producers_or_fixtures_replayed=0,
        source_direction_parameter_precision_cutoff_scan=False,
        full_entry_box_asserted_PSD=False, unseen_far_inverse_energy_upper_supplied=False,
        actual_center_sign_claim=False, actual_full_border_sign_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started,
        inherited_joint_trace_conditions=trace['conditions'],
        inherited_profile_conditions=trace['inherited_profile_conditions'],
        conditions=[
            'All inherited common-form/domain, physical-majorant, positive-reference, source-energy, source/lift and normalization premises remain. Parent-result limitations describe their parent certificates.',
            'The exact W=T R T* and S_F are PSD, including singular cases; no signed-kernel or full interval-entry-box PSD is assumed.',
            'PSD trace monotonicity gives W S_F W <= tr(S_F W) W. Triangle inequalities bound the SAME nominal gradient using saved reference energies and the first coordinate of the SAME conditional mean.',
            'The nominal inverse is W-d e0 e0*, d=-1/p>0. Dropping the negative residual term -r*W r gives an upper for the full nominal curvature.',
            'The saved nominal center, positive nominal denominator, original remaining energy and complete actual-minus-nominal reference/denominator payment are reused without recomputing them.',
            'Every upper is uniform for the same common realization; no component endpoints are asserted jointly attainable.',
            'A negative full upper uses z>=QJ and A>0 to give det(N2)>0 only for the fixed optimistic comparison. It supplies no parameter eligibility, transfer to the actual target, or unseen z upper. A nonnegative upper is inconclusive.',
            'Actual-center/full-border/cofinal positivity, Robin and RH remain unresolved. This is an application of existing PSD trace and norm inequalities, with no generic-method novelty or Lean-proof claim.'])
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
