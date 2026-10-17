"""Compare the same source's negative energy and complete physical cost.

Consume accepted scalar enclosures. The previous norm estimator is unchanged;
no Gram, trace, reference action, inverse, solve or source query is repeated.
Schur elimination, energy completion and quadratic norm expansion are standard.
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
ROOT = 'full-gamma-stationarity-family-norm.json.zlib'
ROOT_SHA = '169e858deb9f06ca677098da032849450074f794627eb0174711717abc5ca1b8'
PRODUCER = 'full_gamma_stationarity_family_norm.py'
FAMILY = 'full-gamma-coherent-family-norm.json.zlib'
CAPACITY = 'full-gamma-negative-channel-capacity.json.zlib'
PROJECTION = 'full-gamma-projection-upper.json.zlib'
PROFILE = 'full-gamma-conditioned-profile.json.zlib'
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'


def joint_budget(d, alpha, schur, beta, n0, gamma, rho_upper, energy_upper,
                 nominal_upper, conditional_upper, inverse_upper, delta):
    """Pointwise SAME-source beta^2 <= (Eplus-alpha)*h is inherited."""
    args = (d, alpha, schur, beta, n0, gamma, rho_upper, energy_upper,
            nominal_upper, conditional_upper, inverse_upper, delta)
    if any(not arb(v).is_finite() for v in args):
        raise ArithmeticError('Nonfinite joint margin or physical-cost data.')
    if (not d > 0 or not schur > 0 or not alpha >= 0 or n0.upper() < 0
            or any(not v >= 0 for v in (rho_upper, energy_upper, nominal_upper,
                                       conditional_upper, inverse_upper, delta))):
        raise ArithmeticError('Missing inherited PSD, energy, denominator or allowance domain.')
    err = (arb(conditional_upper)+arb(inverse_upper)).upper()
    error_cost = (2*arb(nominal_upper).sqrt()*err+err**2).upper()
    coefficient = d-alpha-arb(delta)*(n0.upper()+error_cost)
    signed_gain = -2*arb(delta)*beta*gamma
    energy_cost = (arb(delta)*arb(energy_upper)*arb(rho_upper)).upper()
    omega = coefficient*schur-beta*beta+signed_gain-energy_cost
    # Round the sufficient numerator lower BEFORE choosing its denominator.
    lower = arb(fmpq(int(envelope(omega, BITS)[0]), 1 << BITS))
    denominator = schur.upper() if lower > 0 else schur.lower()
    upper = (-lower/denominator).upper()
    values = dict(negative_pivot_d=d, first_row_energy_alpha=alpha,
        stationarity_schur=schur, second_source_numerator=beta,
        nominal_anchor_norm_squared=n0, signed_nominal_cross=gamma,
        nominal_square_norm_upper=arb(rho_upper), remaining_energy_upper=arb(energy_upper),
        inherited_nominal_norm_squared_upper=arb(nominal_upper),
        inherited_conditional_response_norm_upper=arb(conditional_upper),
        inherited_inverse_reference_error_norm_upper=arb(inverse_upper),
        complete_error_norm_upper=err, complete_error_squared_cost_upper=error_cost,
        complete_budget_coefficient=coefficient, signed_cross_gain=signed_gain,
        conditional_energy_cost_upper=energy_cost, joint_margin_numerator=omega,
        rounded_joint_numerator_lower=lower, complete_joint_family_quadratic_upper=upper)
    if any(not arb(v).is_finite() for v in values.values()):
        raise ArithmeticError('Nonfinite joint comparison result.')
    return values


def consume(out):
    started = time.perf_counter()
    packets, guards = {}, {}
    def read(name, digest):
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted joint-margin input: '+name)
        guards[name] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
        return packets[name]
    root = read(ROOT, ROOT_SHA)
    for name, digest in root['packet_sha256'].items():
        read(name, digest)
    for name, digest in dict(root['inherited_dependency_sha256'],
                             **{PRODUCER:root['source_sha256']}).items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Different accepted joint-margin helper: '+name)
        guards[name] = digest
    family, capacity, projection, profile = (packets[n] for n in (FAMILY, CAPACITY, PROJECTION, PROFILE))
    if (root['schema'] != 'ACTUAL_C27_STATIONARITY_COHERENT_FAMILY_V1'
            or root['source_sha256'] != 'e16031be55706d509027664c7bb0132a7a862e0fc1b1c31caeaec57dd0f28972'
            or root['original_metadata_grid_bits'] != 512 or len(root['remaining_modes']) != 18
            or any(p.get('grid_bits', p.get('inherited_grid_bits')) != BITS
                   or p['remaining_modes'] != root['remaining_modes'] for p in packets.values())
            or any((p['source_integer_sha256'], p['source_bits'], p['source_row'])
                   != (SOURCE_ID, 384, 0) for p in (root, family, capacity, profile))
            or any(p['profile_order'] != root['profile_order'] for p in (family, profile))
            or root['response_order'] != ['negative_column', 'r_reference_response',
                   'bJ_reference_response', 'first_measured_reference_response',
                   'second_measured_reference_response']
            or not profile['conditioning_certified']):
        raise ValueError('Different SAME source, exact frame, original grid or conditional family.')
    rational = lambda p, k: fmpq(*map(int, p['exact_rationals'][k]))
    scalar = lambda p, k: interval_ball(p['scalar_bounds'][k], BITS)
    t_plus, fee, delta = (rational(root, k) for k in
        ('t_plus', 'actual_charge_floor', 'true_baseline_gain_coefficient'))
    if (any(rational(root, k) != rational(family, k) for k in
            ('t_plus', 'actual_charge_floor', 'true_baseline_gain_coefficient'))
            or delta != fmpq(33, 16)*fmpq(profile['exact_far_allowance'])
            or not t_plus > 0 or not fee > 0):
        raise ValueError('Different original endpoint or complete baseline allowance.')
    first_diagonal = scalar(root, 'first_source_diagonal')
    if not first_diagonal > 0:
        raise ArithmeticError('Missing positive first-source diagonal.')
    k1 = interval_ball(projection['matrices']['source_on_lift'][0][0], BITS)
    alpha = k1**2/first_diagonal
    values = joint_budget(-scalar(capacity, 'old_first_pivot'), alpha,
        scalar(root, 'stationarity_schur'), scalar(root, 'stationarity_second_numerator'),
        scalar(root, 'nominal_quadratic_constant'), scalar(root, 'nominal_quadratic_half_linear'),
        max(arb(0), scalar(root, 'nominal_quadratic_square').upper()),
        scalar(root, 'first_row_energy_remainder_upper').upper(),
        scalar(root, 'stationarity_nominal_norm_squared_upper').upper(),
        scalar(root, 'stationarity_conditional_response_norm_upper').upper(),
        scalar(root, 'inherited_inverse_reference_error_norm_upper').upper(), delta)
    exported = {k:envelope(v, BITS) for k,v in values.items()}
    upper = fmpq(int(exported['complete_joint_family_quadratic_upper'][1]), 1 << BITS)
    lower = fmpq(int(exported['rounded_joint_numerator_lower'][0]), 1 << BITS)
    for name, digest in guards.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Accepted input changed during joint comparison: '+name)
    result = dict(schema='ACTUAL_C27_JOINT_MARGIN_PHYSICAL_COST_V1',
        packet_sha256={n:h for n,h in guards.items() if n in packets},
        inherited_dependency_sha256={n:h for n,h in guards.items() if n not in packets},
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_integer_sha256=SOURCE_ID, source_bits=384, source_row=0,
        remaining_modes=root['remaining_modes'], profile_order=root['profile_order'],
        response_order=root['response_order'], grid_bits=BITS, working_bits=ctx.prec,
        original_metadata_grid_bits=512, scalar_bounds=exported,
        exact_rationals={k:[str(v.p),str(v.q)] for k,v in dict(t_plus=t_plus,
            actual_charge_floor=fee, true_baseline_gain_coefficient=delta,
            joint_numerator_lower=lower, complete_joint_family_upper=upper).items()},
        verdict=('same_stronger_comparison_excluded_by_joint_margin'
                 if lower > 0 and upper < 0 else 'inconclusive_joint_margin_budget'),
        new_geometric_contractions_Grams_traces_reference_actions_inverses_solves_source_far_queries=0,
        fixed_vector_selections=0, old_scientific_producers_or_fixtures_replayed=0,
        source_direction_parameter_precision_scan=False, stopped_norm_estimator_changed=False,
        full_entry_box_asserted_PSD=False, actual_retained_positivity_claim=False,
        actual_cofinal_positivity_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started, conditions=[
            'ALL inherited actual/common-source,domain,exact20frame,physical and Euclidean PSD,conditional-source,two-sided reference,complete original512metadata allowance,normalization,endpoint charge and baseline-order premises remain.',
            'For the SAME exact source p=C^-1*k,alpha=k1^2/A,h=b-c^2/A,beta=k2-c*k1/A and p*C*p=alpha+beta^2/h<=Eplus. Thus beta^2<=(Eplus-alpha)*h. h>0 and rho>=0.',
            'For the SAME nominal physical response h0+(beta/h)*h2, n0=h0*H*h0,gamma=h0*H*h2,rho=h2*H*h2. Its norm is bounded by n0+(2beta*gamma+e*rho)/h. Signed beta*gamma is retained,not replaced by its absolute value.',
            'The accepted conditional and inverse-reference norm errors are unchanged. Their complete squared-norm increment is at most Kerr=2sqrt(N0upper)*(econd+einv)+(econd+einv)^2.',
            'The exact old quadratic is -d+alpha+beta^2/h. The stronger endpoint adds at most original_delta times the complete physical norm. Hence complete Q<=-Omega/h with Omega=(d-alpha-delta*(n0+Kerr))*h-beta^2-2delta*beta*gamma-delta*e*rho.',
            'Scalar interval enclosures hold for every SAME source. They need no jointly attainable corners or PSD entry box. Use accepted eupper and rhoupper; round Omega lower before normalization by hupper if positive,otherwise hlower.',
            'Strict complete negativity excludes only SAME dominated stronger endpoint comparison for0<t<=tplus,D>=K>=Dfloor. No actual retained/cofinal positivity,Robin,RH,maximal tolerance or Lean proof is claimed.',
            'One scalar-only application,standard Schur/energy/norm mechanisms,no generic novelty. No stopped-norm refinement or source,direction,parameter,precision or joint-formula retry.'])
    write_packet(out, result)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Use one fresh .zlib output; no joint-budget retry.')
    ctx.prec = 1024
    consume(args.out)


if __name__ == '__main__':
    main()
