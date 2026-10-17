"""Pay the true-base gain on a coherent source-dependent minimizer family.

Reuse the saved border determinant and complete source energies. The new
joint Euclidean trace bounds exact minimizers; no fixed vector is selected,
rounded, inverted or solved, and no accepted scientific producer is replayed.
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
from joint_gamma_boundary import matrix_in, write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
PINS = {
    'full-gamma-conditioned-profile.json.zlib': '2e1c11ef3807445ce6b157d47232bc6086ac96f8a07786e4b03b7355e61c00ae',
    'joint-gamma-6765-rho32-v2-schur.json': '13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd',
    'full-gamma-directional-metric.json.zlib': 'f6bea1b9b43c5ca8c85fc6cd7222038f2a40b81ed1ce2ef873af91ce1ca987a2',
    'full-gamma-finite-good-block-metric.json.zlib': '9ba53b70a0f1cad25b9e92894eff93f1340dc31d1a9bf945c560b6dc6940d072',
    'full-gamma-far-energy-budget.json.zlib': '11ccabd6110436b2b29f1a910c6c36a5fd8479af19d183cf614826ed9ad746c6',
    'full-gamma-conditioned-joint-trace.json.zlib': '6656a9e9ab05dd090168c93a359eeb8930998b89e855755b062f2493634dd0dd',
    'full-gamma-negative-channel-capacity.json.zlib': '1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e',
    'full-gamma-admissible-charge.json.zlib': '9ed4c784f0f38cb430d48c527fac4478a41577069c7431428c08bc27a0e06435',
}
PROFILE, CORE, CERT, METRIC, BUDGET, TRACE, CAPACITY, ADMISSIBLE = PINS


def finite(values):
    if any(not v.is_finite() for v in values):
        raise ArithmeticError('Nonfinite coherent-family geometry or budget.')


def family_geometry(sigma, reference, euclidean, gram, metric_schur, eta):
    """One new joint Euclidean trace; actual PSD forms, not an entry box."""
    n = sigma.nrows()
    if (n < 2 or sigma.ncols() != n
            or (reference.nrows(), reference.ncols()) != (n-1, n-1)
            or (euclidean.nrows(), euclidean.ncols()) != (n+2, n+2)
            or (gram.nrows(), gram.ncols()) != (2, 2)):
        raise ValueError('Incompatible SAME-frame family geometry.')
    finite(list(sigma.entries())+list(reference.entries())+list(euclidean.entries())
           +list(gram.entries())+[arb(metric_schur), arb(eta)])
    pivot, a, c = sigma[0, 0], gram[0, 0], gram[0, 1]
    if not pivot < 0 or not a > 0 or not metric_schur > 0 or not 0 <= eta < 1:
        raise ArithmeticError('Missing negative pivot, positive metric or reference certificate.')
    # The inherited exact two-dimensional Schur formula; no new inverse call.
    inverse_good = arb_mat([[1/a+c*c/(a*a*metric_schur), -c/(a*metric_schur)],
                           [-c/(a*metric_schur), 1/metric_schur]])
    good_trace = sum((euclidean[i, j]*inverse_good[j, i]
                      for i in range(2) for j in range(2)), arb(0))
    transport = arb_mat([[-sigma[0, j+1]/pivot if i == 0 else arb(int(i == j+1))
                          for j in range(n-1)] for i in range(n)])
    vv = arb_mat([[euclidean[i+2, j+2] for j in range(n)] for i in range(n)])
    positive_gram = symmetric(transport.transpose()*vv*transport)
    product = positive_gram*reference
    core_trace = sum((product[i, i] for i in range(n-1)), arb(0))
    negative_norm = euclidean[2, 2]
    if any(v.upper() < 0 for v in (good_trace, core_trace, negative_norm)):
        raise ArithmeticError('The inherited exact PSD geometries have inconsistent enclosures.')
    trace_upper = (max(arb(0), good_trace.upper())
                   +max(arb(0), core_trace.upper())/(1-arb(eta))).upper()
    values = dict(good_euclidean_energy_trace=good_trace,
        reference_positive_euclidean_trace=core_trace,
        actual_joint_euclidean_energy_trace_upper=trace_upper,
        negative_coordinate_euclidean_norm_squared=negative_norm)
    finite(values.values())
    return values


def family_budget(pivot, gram_a, gram_d, measured_v, projection_upper,
                  rho_reference_upper, complete_positive_upper, eta, t_plus,
                  fee, determinant_lower, trace_upper, negative_norm_upper, delta):
    """Normalize the inherited border gap before paying any Euclidean gain."""
    args = (pivot, gram_a, gram_d, measured_v, projection_upper,
            rho_reference_upper, complete_positive_upper, eta, t_plus, fee,
            determinant_lower, trace_upper, negative_norm_upper, delta)
    finite([arb(v) for v in args])
    if (not pivot < 0 or not gram_a > 0 or not gram_d >= 0 or not 0 <= eta < 1
            or not t_plus > 0 or not fee > 0 or not determinant_lower > 0
            or any(not v >= 0 for v in (projection_upper, rho_reference_upper,
                   complete_positive_upper, trace_upper, negative_norm_upper, delta))):
        raise ArithmeticError('Missing coherent-family domain or determinant margin.')
    c11 = (1/arb(t_plus)+gram_a.upper()+arb(rho_reference_upper)/(1-arb(eta))).upper()
    c22 = (arb(fee)+gram_d.upper()+2*measured_v.upper()
           +arb(projection_upper)+arb(complete_positive_upper)).upper()
    if not c11 > 0 or not c22 > 0:
        raise ArithmeticError('Missing positive small-matrix diagonal upper.')
    det_upper = (c11*c22).upper()
    fraction = (arb(determinant_lower)/det_upper).lower()
    if not fraction > 0 or not fraction <= 1:
        raise ArithmeticError('Border/primal determinant relation contradicts its exact domain.')
    d = -pivot
    margin = (d.lower()*fraction).lower()
    positive_energy = (d.upper()*(1-fraction)).upper()
    norm_upper = ((arb(negative_norm_upper).sqrt()
                   +(arb(trace_upper)*positive_energy).sqrt())**2).upper()
    if not margin > 0 or not norm_upper > 0:
        raise ArithmeticError('Missing nonzero primal gap or family normalization.')
    payment = (arb(delta)*norm_upper).upper()
    threshold = (margin/norm_upper).lower()
    upper = (-margin+payment).upper()
    values = dict(positive_source_middle_11_upper=c11,
        positive_source_middle_22_upper=c22, positive_source_middle_determinant_upper=det_upper,
        inherited_border_to_primal_gap_fraction_lower=fraction,
        exact_family_old_negative_margin_lower=margin,
        exact_family_positive_energy_upper=positive_energy,
        exact_family_euclidean_norm_squared_upper=norm_upper,
        certified_true_baseline_coefficient_threshold_lower=threshold,
        complete_true_baseline_gain_upper=payment,
        complete_coherent_family_quadratic_upper=upper)
    finite(values.values())
    return values


def consume(out):
    started = time.perf_counter()
    packets, guarded = {}, {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different saved family input: '+name)
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
        guarded[name] = digest
    profile, core, cert, metric, budget, trace, capacity, admissible = packets.values()
    rational = lambda p, k: fmpq(*map(int, p['exact_rationals'][k]))
    scalar = lambda p, k: interval_ball(p['scalar_bounds'][k], BITS)
    if (any(p.get('grid_bits', p.get('inherited_grid_bits')) != BITS or p['remaining_modes'] != core['remaining_modes']
            for p in packets.values()) or len(core['remaining_modes']) != 18
            or not profile['conditioning_certified'] or not cert['exact_reference_identities_checked']
            or profile['original_metadata_grid_bits'] != 512
            or any((p['source_integer_sha256'], p['source_bits'], p['source_row'])
                   != (SOURCE_ID, 384, 0) for p in (profile, metric, budget, trace, capacity, admissible))):
        raise ValueError('Different actual source, core, original grid or reference certificate.')
    for name in (PROFILE, CORE, CERT, METRIC, BUDGET, CAPACITY):
        if trace['packet_sha256'][name] != PINS[name]:
            raise ValueError('Different accepted complete-energy dependency: '+name)
    if (profile['packet_sha256'][METRIC] != PINS[METRIC]
            or admissible['packet_sha256'][BUDGET] != PINS[BUDGET]
            or admissible['packet_sha256'][CORE] != PINS[CORE]
            or admissible['verdict'] != 'same_old_base_all_valid_charges_excluded'
            or profile['exact_far_allowance'] != trace['exact_far_allowance']
            or metric['positive_modes'] != profile['positive_modes']
            or metric['basis'] != profile['basis']):
        raise ValueError('Different SAME source-family frame, charge or determinant certificate.')
    for name, digest in profile['inherited_dependency_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Inherited helper identity changed: '+name)
        guarded[name] = digest
    eta = rational(budget, 'neumann_eta_upper')
    if eta != fmpq(int(cert['scalar_bounds']['neumann_eta_upper'][1]), 1 << BITS):
        raise ValueError('Different accepted actual/reference error payment.')
    t_plus = rational(admissible, 'comparison_t_plus')
    fee = rational(admissible, 'inherited_actual_common_charge_floor')
    if t_plus != rational(budget, 'comparison_t_plus'):
        raise ValueError('Different optimistic shift.')
    delta = fmpq(33, 16)*fmpq(profile['exact_far_allowance'])
    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], BITS))
    gram = symmetric(matrix_in(metric['joint_inverse_metric'], BITS))
    euclidean = symmetric(matrix_in(profile['joint_euclidean_gram'], BITS))
    reference = arb_mat([[arb(fmpq(v)) for v in row] for row in cert['reference_inverse']])
    geometry = family_geometry(sigma, reference, euclidean, gram,
                               scalar(budget, 'metric_schur'), eta)
    measured = matrix_in(profile['measured_pair'], BITS)
    values = family_budget(sigma[0, 0], gram[0, 0], gram[1, 1], measured[0, 1],
        scalar(budget, 'far_energy_projection_floor').upper(),
        scalar(capacity, 'steering_reference_rRr').upper(),
        scalar(trace, 'actual_positive_core_energy_upper').upper(), eta, t_plus, fee,
        rational(admissible, 'endpoint_comparison_determinant_positive_lower'),
        geometry['actual_joint_euclidean_energy_trace_upper'],
        geometry['negative_coordinate_euclidean_norm_squared'].upper(), delta)
    pairs = {k: envelope(v, BITS) for k, v in dict(**geometry, **values).items()}
    upper = fmpq(int(pairs['complete_coherent_family_quadratic_upper'][1]), 1 << BITS)
    for name, digest in guarded.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Family input changed during application: '+name)
    result = dict(schema='ACTUAL_C27_COHERENT_FAMILY_NORM_V1', packet_sha256=PINS,
        inherited_dependency_sha256=profile['inherited_dependency_sha256'],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_integer_sha256=SOURCE_ID, source_bits=384, source_row=0,
        remaining_modes=core['remaining_modes'], profile_order=profile['profile_order'],
        grid_bits=BITS, working_bits=ctx.prec, original_metadata_grid_bits=512,
        scalar_bounds=pairs, exact_rationals={k: [str(v.p), str(v.q)] for k, v in dict(
            t_plus=t_plus, actual_charge_floor=fee, true_baseline_gain_coefficient=delta,
            complete_family_quadratic_upper=upper).items()},
        verdict=('same_stronger_base_source_family_excluded_by_coherent_norm'
                 if upper < 0 else 'inconclusive_coherent_family_norm_budget'),
        new_joint_euclidean_geometry_trace=1, fixed_vector_selections=0,
        new_reference_actions_inverse_solve_source_far_queries=0,
        old_scientific_producers_or_fixtures_replayed=0,
        source_direction_parameter_precision_cutoff_scan=False,
        full_entry_box_asserted_PSD=False, actual_retained_form_sign_claim=False,
        actual_rh_proof=False, elapsed_seconds=time.perf_counter()-started,
        conditions=[
            'All inherited actual/common-form/domain,SAME source/lift,normalization,PSD physical and Euclidean forms,complete-far allowance,conditional source energies and exact-reference premises remain.',
            'A*Lstar*A=diag(G,Sigma), J=[T,e0], J*Sigma*J=diag(Pplus,-d). The exact positive block is C0=diag(G,Pplus), not the reference inverse. A* (Lbar-Lstar)*A<=delta*A*A.',
            'The restricted good-source metric has22entry G22+2v+QJ. Setting z=QJ in the accepted bound zcrit-QJ<=Uadm is this restricted20frame border; the actual full good metric need not attain its projection lower.',
            'For each SAME actual source choose its exact minimizer p=(Lambda^-1+F*C0^-1*F*)^-1*k, u=-C0^-1*F**p. No common fixed vector or numerical inverse is selected.',
            'detN/detC=1-k**C^-1*k/d. The inherited detNlower is divided by a NEW rounded detCupper with detC<=C11upper*C22upper<=detCupper to obtain a correctly normalized negative primal gap.',
            'The SAME joint Euclidean PSD Gram gives norm(A*D*(u,0))^2<=tr(Hpositive*C0^-1)*(u*C0*u). Trace monotonicity with R/(1-eta) bounds only the positive inverse. Triangle then bounds norm(A*D*(u,1)); cross terms are paid, not assumed zero.',
            'A strict negative complete upper excludes the SAME stronger lower comparison for0<t<=tplus,D>=K>=Dlo. A nonnegative upper is inconclusive and triggers no direction,parameter,precision,source or norm-refinement retry.',
            'The positive coefficient threshold bounds only this sufficient family certificate. It is not a maximal tolerable perturbation, actual spectral gap, stronger-form sign, actual retained/cofinal sign,Robin,RH or Lean proof. Existing determinant,Schur,trace,energy and norm methods only.'])
    write_packet(out, result)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Use one fresh .zlib output; no family-budget retry.')
    ctx.prec = 1024
    consume(args.out)


if __name__ == '__main__':
    main()
