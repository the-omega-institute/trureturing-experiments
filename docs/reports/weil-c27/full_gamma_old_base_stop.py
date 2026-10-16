"""Test the true old-core, fixed-charge recipe before filling its good metric.

Reuse the accepted17-coordinate trial and its unclipped lower value.
Only its first gradient and one mixed scalar contraction are missing.
The comparison t_plus need not itself be an admitted coefficient.
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
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
MODES = [0, 1, 2, 3, 4, 5, 6, 9, 10, 12, 14, 15, 18, 19, 24, 29, 33, 39]
PINS = {
    'full-gamma-actual-joint-lower.json.zlib':
        'f294428dc9685e04fca777a5dc71ac0374fc17a4c69246422e8b4325351da5a9',
    'joint-gamma-6765-rho32-v2-schur.json':
        '13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd',
    'full-gamma-projection-upper.json.zlib':
        'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3',
    'full-gamma-actual-channel-projection.json.zlib':
        '63e3c2e291940ca3549fa4a1303b55b4d92c08d594a50daf7ca2132087ad2745',
    'full-gamma-surplus-budget.json':
        '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f',
    'full-gamma-band-surplus-budget.json.zlib':
        '04fb6eb533b8448ade52ad4032faaaf079ad6603f9998d272a36828c57be44b1',
    'joint-gamma-6765-rho32-v2-rank-floor.json':
        'f4ceb4b526de202b540089f7532ac3aee71e7013df98a7dfe0ab3975345cdf89',
}


def transport(a, c, beta, rho, x, q0_lower, t_plus, charge):
    """Real fixed-trial identity; no inverse, solve or trial quadratic."""
    if len(c) != len(x) or len(beta) != len(x)+1 or len(rho) != len(beta):
        raise ValueError('Incompatible real core/trial row dimensions.')
    if any(not v.is_finite() for v in [a, q0_lower, charge, t_plus]+c+beta+rho+x):
        raise ArithmeticError('Nonfinite input; stop without retry.')
    if not a < 0 or not t_plus > 0 or not charge > 0:
        raise ArithmeticError('Require a negative old pivot and positive comparison/charge.')
    updated_pivot = a+t_plus*rho[0]**2
    eta = 1+t_plus*rho[0]**2/a
    if not updated_pivot < 0 or not eta > 0:
        raise ArithmeticError('The updated negative-pivot/eta guard is inconclusive.')
    first_gradient = (beta[0]-sum((u*v for u, v in zip(c, x)), arb(0)))/a
    mixed = rho[0]*first_gradient+sum((u*v for u, v in zip(rho[1:], x)), arb(0))
    loss_upper = (t_plus*mixed**2/eta).upper()
    trial_lower = (q0_lower-loss_upper).lower()
    margin_lower = (charge+trial_lower).lower()
    result = dict(old_first_pivot=a, updated_first_pivot=updated_pivot,
                  positive_eta=eta, first_trial_gradient=first_gradient,
                  mixed_source_trial=mixed, cached_trial_lower=q0_lower,
                  transport_loss_upper=loss_upper, updated_trial_lower=trial_lower,
                  fixed_charge=charge, comparison_coefficient=t_plus,
                  stopping_margin_lower=margin_lower)
    if any(not v.is_finite() for v in result.values()):
        raise ArithmeticError('Nonfinite scalar output; stop without retry.')
    return result


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            invoked.add(Path(filename).resolve())
    guards = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    packets = []
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted input: '+name)
        packets.append(json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw))
    joint, core, projection, actual, budget, band, rank = packets
    if [p['schema'] for p in packets] != [
            'ACTUAL_C27_FULL_GAMMA_ACTUAL_JOINT_LOWER_V1',
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1',
            'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1',
            'ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
            'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1',
            'ACTUAL_C27_BAND_SURPLUS_BUDGET_V1',
            'ACTUAL_C27_SMALL_SCHUR_RANK_FLOOR_V1']:
        raise ValueError('Different cached packet roles.')
    core_pin = PINS['joint-gamma-6765-rho32-v2-schur.json']
    actual_pin = PINS['full-gamma-actual-channel-projection.json.zlib']
    projection_pin = PINS['full-gamma-projection-upper.json.zlib']
    if (joint['packet_sha256']['core'] != core_pin
            or joint['packet_sha256']['actual'] != actual_pin
            or projection['packet_sha256']['core'] != core_pin
            or rank['schur_result_sha256'] != core_pin
            or budget['packet_sha256']['full-gamma-projection-upper.json.zlib'] != projection_pin
            or band['packet_sha256']['full-gamma-projection-upper.json.zlib'] != projection_pin
            or band['packet_sha256']['full-gamma-surplus-budget.json'] != PINS['full-gamma-surplus-budget.json']):
        raise ValueError('Different true old core, trial source or band inputs.')
    for role in ('shared', 'projection', 'common', 'defect'):
        if joint['packet_sha256'][role] != actual['packet_sha256'][role]:
            raise ValueError('Different cached finite source: '+role)
    for role in ('projection', 'defect'):
        if projection['packet_sha256'][role] != actual['packet_sha256'][role]:
            raise ValueError('Different unweighted source/lift: '+role)
    target, pivots = core['target_packet_sha256'], core['pivot_result_sha256']
    if (target != '7546a915368b4a40e0e41ea0ea693403c239598c7a6e96b7b8b07e67e22e98c6'
            or pivots != 'd441b9f9fc49338a121430a501043b76afcf7dffa9ec09c8d672aa9a2358b3f0'
            or rank['target_packet_sha256'] != target or rank['pivot_result_sha256'] != pivots
            or projection['packet_sha256']['base'] != target
            or actual['packet_sha256']['base'] != target):
        raise ValueError('Different old retained target or positive-block certificate.')
    expected = core['input_sha256']
    if (any(p['input_sha256'] != expected or p['remaining_modes'] != MODES
            or p['grid_bits'] != BITS for p in (joint, core, projection, actual))
            or rank['remaining_modes'] != MODES
            or any((p['N'], p['M']) != (1280, 32768) for p in (core, projection, actual, budget))
            or core['last_mode'] != 2048 or actual['first_mode'] != 2049
            or core['basis'] != '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'):
        raise ValueError('Different source, real basis or ordered18-core.')
    source_id = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
    if (any((p['source_bits'], p['source_row'], p['source_integer_sha256']) != (384, 0, source_id)
            for p in (budget, band))
            or projection['exact_source_bits'][0] != 384
            or projection['independent_source_indices'][0] != 0
            or projection['exact_source_integer_sha256'][0] != source_id
            or not band['complete_budget_below_conservative_capacity_target']
            or band['verdict'] != 'conditional-band-budget-below-projection-capacity-target'):
        raise ValueError('Different original source or supplied eligible band.')
    minors = rank['leading_minors']
    if (rank['verdict'] != 'full-complement-inertia-counted' or rank['negative_count_lower'] != 1
            or rank['reused_positive_coordinate_count'] != 1263 or len(minors) != 18
            or any((v['order'], v['determinant_sign'], v['pivot_sign'], v['negative_count_through_order'])
                   != (j+1, -1, -1 if j == 0 else 1, 1) for j, v in enumerate(minors))):
        raise ValueError('Require the saved true-core inertia certificate; no determinant replay.')
    producers = {}
    for packet, name in zip(packets, ('full_gamma_actual_joint_lower.py',
            'structured_residual_schur_finish.py', 'full_gamma_projection_upper.py',
            'full_gamma_actual_channel_projection.py', 'full_gamma_surplus_budget.py',
            'full_gamma_band_surplus_budget.py', 'small_schur_rank_floor.py')):
        producer = directory/name
        digest = hashlib.sha256(producer.read_bytes()).hexdigest()
        if digest != packet['source_sha256']:
            raise ValueError('Different cached producer source: '+producer.name)
        producers[producer] = digest
    clipped_lower = joint['scalar_bounds']['new_positive_energy_lower']
    if int(clipped_lower[0]) <= 0:
        raise ArithmeticError('Cached max(0,trial) clipping must be inactive for this exact trial.')
    sigma = core['reduced_schur_intervals']
    beta_raw = actual['actual_finite_row_on_lift']
    rho_raw = projection['matrices']['source_on_lift']
    profile = joint['new_reference_profile']
    if (len(sigma) != 18 or any(len(row) != 18 for row in sigma)
            or len(beta_raw) != 1 or len(beta_raw[0]) != 18
            or len(rho_raw) != 2 or any(len(row) != 18 for row in rho_raw)
            or len(profile) != 17 or any(len(row) != 1 for row in profile)):
        raise ValueError('Incomplete accepted first row, source row or exact trial.')
    a = interval_ball(sigma[0][0], BITS)
    c = [interval_ball(pair, BITS) for pair in sigma[0][1:]]
    if any(not c[j-1].overlaps(interval_ball(sigma[j][0], BITS)) for j in range(1, 18)):
        raise ArithmeticError('Incompatible true-core first-row symmetric enclosures.')
    beta = [interval_ball(pair, BITS) for pair in beta_raw[0]]
    rho = [interval_ball(pair, BITS) for pair in rho_raw[0]]
    x = [arb(fmpq(row[0])) for row in profile]
    q0 = fmpq(int(joint['scalar_bounds']['actual_finite_signed_energy_lower'][0]), 1 << BITS)
    charge = fmpq(int(actual['scalar_bounds']['refined_common_charge_upper'][1]), 1 << BITS)
    hlo = fmpq(*map(int, budget['exact_rationals']['h_lower']))
    hup = fmpq(*map(int, budget['exact_rationals']['h_upper']))
    cap = fmpq(*map(int, band['exact_rationals']['supplied_parameter_cap']))
    if not 0 < hlo <= hup or not cap > 0:
        raise ArithmeticError('Require positive ordered h endpoints and supplied cap.')
    t_plus = cap/hlo
    scalars = transport(a, c, beta, rho, x, arb(q0), arb(t_plus), arb(charge))
    exported = {name: envelope(value, BITS) for name, value in scalars.items()}
    if any(pair is None or int(pair[0]) > int(pair[1]) for pair in exported.values()):
        raise ArithmeticError('Invalid outward scalar export; stop without retry.')
    trigger = scalars['stopping_margin_lower'] > 0
    result = dict(schema='ACTUAL_C27_TRUE_OLD_BASE_STOP_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={p.name: h for p, h in guards.items()},
        cached_producer_source_sha256={p.name: h for p, h in producers.items()},
        input_sha256=expected, target_packet_sha256=target, pivot_result_sha256=pivots,
        source_integer_sha256=source_id, source_bits=384, source_row=0,
        N=1280, M=32768, remaining_modes=MODES, grid_bits=BITS, working_bits=1024,
        exact_rationals={name: [str(v.p), str(v.q)] for name, v in dict(
            h_lower=hlo, h_upper=hup, supplied_parameter_cap=cap,
            upper_comparison_coefficient=t_plus, fixed_charge=charge, cached_trial_lower=q0).items()},
        scalar_bounds=exported, cached_trial_clipping_inactive=True,
        comparison_coefficient_eligibility_claim=False,
        verdict=('fixed-charge-zero-response-old-base-recipe-excluded-through-supplied-band'
                 if trigger else 'inconclusive-true-old-base-stop'),
        new_first_gradient_terms=17, new_mixed_contraction_terms=18,
        old_producers_replayed=0, new_inverse_or_solve_calls=0,
        new_trial_or_quadratic_forms=0, new_source_or_gram_actions=0,
        parameter_source_precision_scan=False, elapsed_seconds=time.monotonic()-started,
        conditions=joint['conditions']+band['conditions']+rank['conditions']+[
            'The cached finite signed lower is for the SAME exact17-trial: the saved positive lower is strictly positive, so max(0,variational.lower()) did not clip.',
            'Sigma is the TRUE old Schur core with the saved one-negative inertia; SigmaMax and other compression uppers are not substituted.',
            't_plus=tau_cap/h_lower is only an upper comparison, dominating actual tau/h for every supplied 0<=tau<=tau_cap; its own eligibility is not asserted.',
            'For zero projected far response, a strict trigger makes the center border positive definite for any PSD joint good-block metric; the corresponding auxiliary target has a negative direction for 0<=t<=t_plus,D>=fixed_charge.',
            'The robust obstruction applies only to an explicitly enlarged energy/Gram upper-bound relaxation containing zero projected far response. No such completion is identified with the actual fixed source.'],
        scope='Fixed-charge, true-old-base, zero-projected-far-response auxiliary recipe only. '
            'Standard cached variational transport and inertia application, not a new general theorem. '
            'Actual source far response, smaller source-dependent D, stronger true baseline, '
            'retained/cofinal positivity, Robin tail, RH and Lean remain unresolved.')
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Accepted input changed: '+name)
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Invoked source/helper changed: '+path.name)
    for path, digest in producers.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Cached producer source changed: '+path.name)
    with out.open('xb') as stream:
        stream.write(zlib.compress(json.dumps(result, separators=(',', ':')).encode(), 9))
    for name, value in scalars.items():
        print(name, value.str(20), flush=True)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Use a fresh .zlib output; accepted archives are never overwritten.')
    ctx.prec = 1024
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
