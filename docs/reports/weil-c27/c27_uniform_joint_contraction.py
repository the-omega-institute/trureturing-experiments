"""Consume uniform actual slack through one whole-inverse comparison.

Use only saved rational endpoints and all exterior diagonal weights. No
actual direction, trial quadratic, residual, Gram, inverse or solve runs.
A negative value rejects this sufficient target, not the actual Weil form.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time
import zlib

HERE = Path(__file__).resolve().parent
PINS = {'propagated-prime-row.json.zlib': '935b2c0707c3187fa652730b1e40d5937a4e071158b7cb9f1f4a6a78f8459ce6', 'inverse-residual-pilot.json.zlib': '6db5093bf78fe1dcf43cbbe1b39ef2b1faf4174cf7cf33fe796c0d10ad4d231b', 'joint-gamma-6765-direction.json': 'ccdb3360923f4fd07086871c5d02512a71301143376db8fd54aa23a3797152cc', 'full-gamma-dual-32768-finite.json.zlib': '4085fa7048bafc987c75b2340b56d98d556534ae74d304e61469d4e6e52af1b5'}
ORIGINAL = 'ea4a761f76b59a11b2d04faf80531901a5164f4145c2424a4e79889ae3f72e22'
BASE = '7546a915368b4a40e0e41ea0ea693403c239598c7a6e96b7b8b07e67e22e98c6'
BASIS = '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'


def pair(q):
    return [str(q.numerator), str(q.denominator)]


def endpoints(value, bits):
    if (type(bits) is not int or bits <= 0 or not isinstance(value, list)
            or len(value) != 2 or any(type(x) is not str for x in value)):
        raise ValueError('Require ordered exact endpoint strings and a positive grid.')
    lo, hi = (F(int(x), 1 << bits) for x in value)
    if lo > hi:
        raise ValueError('Reversed saved interval.')
    return lo, hi


def exterior_bound(pilot):
    """Bound Lambda on ALL exterior modes, including selected E and the far band."""
    bits = pilot['middle_weight_bits']
    nums = pilot['middle_weight_numerators']
    if (type(bits) is not int or bits <= 0 or len(nums) != 32768-1280
            or any(type(x) is not str for x in nums)):
        raise ValueError('Require every original1281..32768 weight.')
    middle = [F(int(x), 1 << bits) for x in nums]
    far = F(int(pilot['far_weight_numerator']), int(pilot['far_weight_denominator']))
    if far <= 0 or any(x <= 0 for x in middle):
        raise ValueError('Nonpositive exterior weight; no reciprocal allowed.')
    minimum = min(middle+[far])
    location = 1281+middle.index(minimum) if minimum in middle else 'far>=32769'
    lam = 1/minimum
    # G2, pi<22/7, exp(2)>7, sqrt(27)>5 and log(27)<4.
    # 2*pi/(1-exp(-2)) < 22/3;
    # 2*d/((sqrt27-1)*(1-1/729)) < 729/364.
    # Complete even pole norm <=4*sinh(log27/2)=52/sqrt27<52/5.
    gamma, pole = F(22, 3)+F(729, 364), F(52, 5)
    return dict(minimum_weight=minimum, maximum_lambda=lam,
                maximum_lambda_location=location, gamma_norm_upper=gamma,
                pole_norm_upper=pole, comparison_norm_upper=lam+gamma+pole)


def contract_cached(q, old, gap, norm):
    """An exact necessary test, never a full-matrix positivity test."""
    qlo, qhi = q
    olo, ohi = old
    if qlo > qhi or olo > ohi or qlo <= 0 or ohi >= 0 or gap <= 0 or norm <= 0:
        raise ValueError('Need positive Q, negative old value and positive gap/norm.')
    c = norm/(norm+gap)
    value = ((norm*olo+gap*qlo)/(norm+gap), (norm*ohi+gap*qhi)/(norm+gap))
    ceiling = (gap*qlo/(-olo), gap*qhi/(-ohi))
    verdict = ('fixed_uniform_contraction_insufficient' if value[1] < 0 else
               'necessary_direction_positive_only' if value[0] > 0 else 'inconclusive')
    return c, value, ceiling, verdict


def consume(directory, out):
    started = time.perf_counter()
    packets = []
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted input: '+name)
        packets.append(json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw))
    row, pilot, direction, finite = packets
    if [p['schema'] for p in packets] != [
            'ACTUAL_C27_CANONICAL_PROPAGATED_PRIME_ROW_V1',
            'ACTUAL_C27_INVERSE_RESIDUAL_PILOT_V1',
            'ACTUAL_C27_JOINT_GAMMA_DIRECTION_V1',
            'ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1']:
        raise ValueError('Different accepted packet roles.')
    meta = finite['metadata']
    if (meta['direction_sha256'] != PINS['joint-gamma-6765-direction.json']
            or meta['target_sha256'] != ORIGINAL or direction['original_target_sha256'] != ORIGINAL
            or meta['input_sha256']['inverse-residual-pilot.json.zlib'] != PINS['inverse-residual-pilot.json.zlib']
            or direction['updated_target_sha256'] != BASE
            or finite['boundary_sha256'] != direction['finite_action_sha256']
            or finite['boundary_sha256'] != '94aec57dacacad6bae1698a7e643514e3dace976fddbcc0008b2b7fd5aab447c'
            or any((p['N'], p['M']) != (1280, 32768) for p in [pilot, direction, meta])
            or meta['last_mode'] != 2048 or direction['selected_remaining_mode'] != 0
            or any(p['basis'] != BASIS for p in [pilot, direction, meta])
            or direction['direction_dyadic_bits'] != 256
            or len(direction['direction_dyadic_integers']) != 1281
            or direction['grid_bits'] != 640 or finite['grid_bits'] != 640
            or finite['dual_bits'] != 384
            or row['canonical_weight'] != 'q1=C_p*(1+(29/50)*(x/L)^2)'
            or row['propagation_applications'] != 1
            or row['verdict'] != 'actual_slack_and_new_eligible_coefficient'
            or row['actual_retained_or_cofinal_positivity_claim'] is not False):
        raise ValueError('Different whole direction, original trial, exterior, basis or surplus.')
    bound = exterior_bound(pilot)
    gap = F(*map(int, row['exact_rationals']['actual_surplus_floor_lower']))
    if gap <= F(9, 25):
        raise ValueError('Missing accepted positive actual uniform surplus.')
    q = endpoints(finite['Q'], finite['grid_bits'])
    old = endpoints(direction['common_young_target_quadratic'], direction['grid_bits'])
    c, value, ceiling, verdict = contract_cached(q, old, gap, bound['comparison_norm_upper'])
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Input changed during consumption: '+name)
    rational = {k: pair(v) for k, v in bound.items() if isinstance(v, F)}
    rational.update(uniform_gap_lower=pair(gap), contraction_factor=pair(c),
                    contraction_gain_fraction=pair(1-c))
    intervals = {k: [pair(lo), pair(hi)] for k, (lo, hi) in dict(
        initial_trial_Q=q, old_young_target=old, contracted_target=value,
        necessary_comparison_norm_ceiling=ceiling).items()}
    result = dict(schema='ACTUAL_C27_UNIFORM_JOINT_CONTRACTION_V1',
        packet_sha256=PINS, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        N=1280, M=32768, basis=BASIS, last_mode=2048,
        original_target_sha256=ORIGINAL, base_target_sha256=BASE,
        direction_packet_sha256=PINS['joint-gamma-6765-direction.json'],
        direction_bits=256, retained_direction_coordinates=1281,
        source_Q_is_scalar_only=True, exact_rationals=rational, exact_intervals=intervals,
        maximum_lambda_location=bound['maximum_lambda_location'],
        verdict=verdict, fixed_comparison_applications=1,
        old_scientific_producers_or_fixtures_replayed=0,
        new_direction_trial_residual_Gram_inverse_solve_far_or_full_target_actions=0,
        norm_weight_direction_parameter_precision_scan=False,
        actual_retained_or_cofinal_positivity_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started,
        conditions=[
            'Same actual common-form/domain comparison T=Lambda+B_GP>0 and H>=T+gI; all Gamma rates and positive even pole retained.',
            'G2 bounds the COMPLETE Gamma boundary; compression/restriction do not increase norm. Pole norm<=4*sinh(L). The rational replacements bound from above, without claiming a sharp norm.',
            'Lambda=W^-1 uses ALL archived exterior weights, including selected E and constant far D, not a J-only minimum.',
            'The original19-channel Young comparison Wbar dominates T^-1 on the residual range, and the SAME initial-trial target satisfies Lstar<=Q-R*Wbar R.',
            'Inverse order gives H^-1<=c*T^-1<=c*Wbar and the original residual identity gives S_actual>=c*Lstar+(1-c)*Q on ALL retained coefficients.',
            'Only Q(v) and Lstar(v) of the SAME complete exact256-bit direction are consumed. A negative upper rejects this fixed sufficient comparison; no fullQ diagonal is recovered and no actual negative direction follows.',
            'Inherited exact trial, normalization, source and sound-cache premises remain. No generic novelty, Lean verification, actual retained/cofinal positivity, Robin or RH conclusion.',
        ])
    with out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(verdict, 'seconds', result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.json':
        parser.error('Require a fresh JSON output; no overwrite or retry.')
    consume(HERE, args.out)
