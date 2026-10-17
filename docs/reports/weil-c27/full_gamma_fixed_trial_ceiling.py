"""Screen refinement of the one saved common-source covector.

Reuse exact convex-endpoint logic; no physical-source supplier is called.
The bound limits an unchanged lower-certificate recipe, not actual energy.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time
import zlib

from full_gamma_response_ceiling import optimistic_ceiling, read_rational, encode

HERE = Path(__file__).resolve().parent
INPUT = 'full-gamma-common-trial-residual.json.zlib'
INPUT_SHA = 'c0aa04f3fbf6d523d0846003a8d2b8385a3d6d3fc891a882d14953a9e88f27c7'
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'


def fixed_trial_capacity(cross_lo, cross_hi, norm_lo, norm_hi):
    """Cap dist(0,J)^2/U for J within I and norm_lo<=U<=norm_hi."""
    if any(type(x) is not Fraction for x in (cross_lo, cross_hi, norm_lo, norm_hi)):
        raise TypeError('Exact Fraction interval endpoints required.')
    if cross_lo > cross_hi or norm_lo > norm_hi:
        raise ValueError('Ordered saved intervals required.')
    if norm_lo <= 0:
        raise ValueError('Ceiling inconclusive without a positive dual-norm lower.')
    return optimistic_ceiling(1/norm_lo, Fraction(0), Fraction(0), cross_lo, cross_hi)[2]


def consume(out):
    started = time.perf_counter()
    source = Path(__file__).resolve()
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    raw = (HERE/INPUT).read_bytes()
    if hashlib.sha256(raw).hexdigest() != INPUT_SHA:
        raise ValueError('Different accepted common-trial input.')
    data = json.loads(zlib.decompress(raw))
    if (data['schema'] != 'ACTUAL_C27_COMMON_TRIAL_RESIDUAL_V1'
            or (data['source_bits'], data['source_row'], data['grid_bits'],
                data['working_bits'], data['source_integer_sha256']) !=
               (384, 0, 640, 1024, SOURCE_ID)
            or any(data[k] for k in ('actual_paid_energy_upper_claim',
                'actual_center_sign_claim', 'actual_full_border_sign_claim', 'actual_rh_proof'))):
        raise ValueError('Different source, form scope or saved claim boundary.')
    intervals = {}
    for key in ('projected_residual_cross', 'projected_dual_norm'):
        pair = data['scalar_bounds'][key]
        if (not isinstance(pair, list) or len(pair) != 2
                or any(type(x) is not str for x in pair)):
            raise ValueError('Expected two saved dyadic endpoints: '+key)
        intervals[key] = tuple(Fraction(int(x), 1 << data['grid_bits']) for x in pair)
    cross_lo, cross_hi = intervals['projected_residual_cross']
    norm_lo, norm_hi = intervals['projected_dual_norm']
    capacity = fixed_trial_capacity(cross_lo, cross_hi, norm_lo, norm_hi)
    rat = lambda key: read_rational(data['exact_rationals'], key)
    if not 0 <= rat('new_residual_lower') <= capacity:
        raise ArithmeticError('Fixed-trial ceiling contradicts saved same-direction residual lower.')
    recipes = {}
    for recipe in ('fixed', 'energy_endpoint'):
        same_F, previous = rat(recipe+'_same_F_lower'), rat(recipe+'_retained_lower')
        ceiling = max(previous, same_F+capacity)
        recipes[recipe] = dict(
            same_F_lower=encode(same_F), current_retained_lower=encode(previous),
            fixed_trial_retained_certificate_ceiling=encode(ceiling),
            fixed_trial_refinement_cannot_certify_positive=ceiling < 0)
        if same_F < 0:
            recipes[recipe].update(
                cross_squared_needed_for_positive=encode(-same_F*norm_lo),
                available_capacity_over_needed_residual=encode(capacity/(-same_F)))
    result = dict(schema='C27_FIXED_TRIAL_CERTIFICATE_CEILING_V1',
        input_packet=INPUT, input_sha256=INPUT_SHA, source_sha256=source_sha,
        reused_endpoint_helper='full_gamma_response_ceiling.optimistic_ceiling',
        source_integer_sha256=SOURCE_ID, source_bits=384, source_row=0,
        grid_bits=640, working_bits='exact rational',
        scalar_intervals={k:[encode(x) for x in v] for k,v in intervals.items()},
        fixed_trial_residual_certificate_ceiling=encode(capacity), recipes=recipes,
        both_fixed_trial_refinements_excluded=all(
            v['fixed_trial_refinement_cannot_certify_positive'] for v in recipes.values()),
        new_saved_fixed_trial_capacity_applications=1, new_recipe_certificate_ceilings=2,
        old_producers_or_fixtures_replayed=0,
        new_physical_source_reference_inverse_solve_or_tail_calls=0,
        source_direction_parameter_precision_cutoff_scan=False,
        actual_paid_energy_upper_claim=False, actual_center_sign_claim=False,
        actual_full_border_sign_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started, conditions=[
            'All original source, Q, physical-error, plane and covector bindings remain. This is an application to the saved common-trial intervals, not a new source measurement.',
            'For any nonempty closed J inside the saved cross interval, and any positive dual-norm upper U obtained from a nonempty refined norm interval inside the saved interval, dist(0,J)^2/U is at most max(abs(Ilo),abs(Ihi))^2/Dlo. Dlo>0 is required only for this ceiling, not for the original residual lower.',
            'Each SAME-F lower and its previous retained lower are held fixed. Maximum retention, never addition of lower bounds on the same residual, gives the reported recipe ceiling.',
            'A negative ceiling rules out positivity only by refining this fixed covector within these intervals with the SAME-F certificate unchanged. It does not upper-bound total residual, actual paid energy or actual source sign.',
            'Joint F/residual estimates, a new covector or a stronger SAME-F lower are outside this diagnostic. RH, Robin, actual-center/full-border/good-block/cofinal positivity remain unresolved; no generic novelty or Lean-proof claim.',
        ])
    with Path(out).open('xb') as stream:
        stream.write(zlib.compress(json.dumps(result,sort_keys=True,separators=(',',':')).encode()))
    print('saved fixed-trial ceiling supplied; seconds',result['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Require a fresh .zlib output; never overwrite accepted data.')
    consume(args.out)
