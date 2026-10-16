"""Bound the same archived action's surplus budget using cached energies.

Only exact rational scalar arithmetic is new. No source/action/tail,
inverse, Gram, projection, parameter search or precision retry is run.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time
import zlib


HERE = Path(__file__).resolve().parent
PINS = {
    'full-gamma-common-column.json.zlib':
        '314e341189f72516aa9f251c6b873587f1bfae54d902c35f6c87c9b03130dfe7',
    'inverse-residual-pilot.json.zlib':
        '6db5093bf78fe1dcf43cbbe1b39ef2b1faf4174cf7cf33fe796c0d10ad4d231b',
    'full-gamma-dual-32768-finite.json.zlib':
        '4085fa7048bafc987c75b2340b56d98d556534ae74d304e61469d4e6e52af1b5',
    'full-gamma-dual-32768-defect.json.zlib':
        '6f320576338562a4e7338906ca66b25ab49e4a1d4460e87e26165cd9579a1154',
    'full-gamma-primal-tail.json':
        'fad634fb25d8d20e83a15d8277569978646246c2a560d0dc0f7a854b5087455a',
    'full-gamma-projection-upper.json.zlib':
        'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3',
}


def rational(value):
    return [str(value.numerator), str(value.denominator)]


def interval(pair, bits):
    lo, hi = (Fraction(int(x), 1 << bits) for x in pair)
    if lo > hi:
        raise ValueError('Reversed accepted scalar interval.')
    return lo, hi


def consume(directory, out):
    started = time.monotonic()
    packets = {}
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted input: '+name)
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    common, pilot, finite, defect, primal, projection = packets.values()
    if [p['schema'] for p in packets.values()] != [
            'ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1',
            'ACTUAL_C27_INVERSE_RESIDUAL_PILOT_V1',
            'ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1',
            'ACTUAL_C27_FINITE_DUAL_DEFECT_V1',
            'ACTUAL_C27_FULL_GAMMA_PRIMAL_TAIL_V1',
            'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1']:
        raise ValueError('Different accepted packet roles.')
    metadata = finite['metadata']
    if any(p['metadata'] != metadata for p in (common, defect, primal)):
        raise ValueError('Different common source/action metadata.')
    if ((metadata['N'], metadata['M'], metadata['last_mode']) != (1280, 32768, 2048)
            or any((p['dual_first'], p['dual_last']) != (2049, 32768)
                   for p in (finite, defect))
            or (pilot['N'], pilot['M']) != (1280, 32768)
            or any(p['grid_bits'] != 640 for p in (common, defect, primal, projection))):
        raise ValueError('Different support or accepted scalar grid.')
    for packet in (common, primal):
        if (packet['finite_action_sha256'] != PINS['full-gamma-dual-32768-finite.json.zlib']
                or packet['finite_defect_sha256'] != PINS['full-gamma-dual-32768-defect.json.zlib']
                or packet['pilot_sha256'] != PINS['inverse-residual-pilot.json.zlib']
                or packet['direction_sha256'] != metadata['direction_sha256']):
            raise ValueError('Different full action or original far weight.')
    if (defect['finite_action_sha256'] != PINS['full-gamma-dual-32768-finite.json.zlib']
            or common['primal_tail_sha256'] != PINS['full-gamma-primal-tail.json']
            or metadata['input_sha256']['inverse-residual-pilot.json.zlib'] !=
            PINS['inverse-residual-pilot.json.zlib']
            or projection['packet_sha256']['finite'] != PINS['full-gamma-dual-32768-finite.json.zlib']
            or projection['packet_sha256']['defect'] != PINS['full-gamma-dual-32768-defect.json.zlib']
            or projection['input_sha256'] != metadata['input_sha256']
            or (projection['N'], projection['M'], projection['finite_first'], projection['finite_last']) !=
            (1280, 32768, 2049, 32768)):
        raise ValueError('Different cached energy/h provenance.')
    integers = finite['dual_integers']
    source_hash = hashlib.sha256(('\n'.join(str(int(x)) for x in integers)+'\n').encode()).hexdigest()
    if (len(integers) != 30720 or finite['dual_bits'] != 384
            or projection['exact_source_bits'][0] != 384
            or projection['exact_source_integer_sha256'][0] != source_hash
            or source_hash != '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'):
        raise ValueError('Different original unscaled w0.')

    # W_J is indexed from 1281; only its 2049:32768 slice acts on u_J.
    weights = list(map(int, pilot['middle_weight_numerators']))
    if len(weights) != 31488 or pilot['middle_weight_bits'] != 256:
        raise ValueError('Incomplete original middle weight array.')
    wmin_numerator = min(weights[768:])
    wmin = Fraction(wmin_numerator, 1 << 256)
    dfar = Fraction(int(pilot['far_weight_numerator']), int(pilot['far_weight_denominator']))
    ej_lo, ej_up = interval(common['bounds']['middle_action_energy'], 640)
    ef_lo, ef_up = interval(common['bounds']['far_action_energy_upper'], 640)
    h_lo, h_up = interval(projection['matrices']['source_gram'][0][0], 640)
    if min(wmin, dfar, h_lo) <= 0 or ej_lo < 0 or ef_lo < 0:
        raise ValueError('Require certified positive weights/h and nonnegative energies.')
    # The archived far energy is Dfar * ||u_far||^2 with CONSTANT Dfar.
    middle_norm_upper, far_norm_upper = ej_up/wmin, ef_up/dfar
    norm_upper = middle_norm_upper+far_norm_upper
    if norm_upper <= 0:
        raise ValueError('The nonzero common action must have a positive norm upper.')
    sigma = Fraction(1, 632)
    p_upper = norm_upper/sigma
    t_cap = h_lo/(h_lo+p_upper)
    exponent = max(1, t_cap.denominator.bit_length()-t_cap.numerator.bit_length()+2)
    t_example = Fraction(1, 1 << exponent)
    if not 0 < t_example <= t_cap < 1 or p_upper > h_lo*(1-t_example)/t_example:
        raise ArithmeticError('The exact example fails the directional budget.')
    result = dict(
        schema='ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1', packet_sha256=PINS,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        N=1280, M=32768, finite_first=2049, finite_last=32768,
        source_row=0, source_bits=384, source_integer_sha256=source_hash,
        exact_rationals={name: rational(value) for name, value in dict(
            minimum_middle_weight=wmin, constant_far_weight=dfar,
            middle_action_energy_upper=ej_up, far_action_energy_upper=ef_up,
            middle_unweighted_norm_squared_upper=middle_norm_upper,
            far_unweighted_norm_squared_upper=far_norm_upper,
            full_unweighted_norm_squared_upper=norm_upper,
            h_lower=h_lo, h_upper=h_up, global_surplus_floor=sigma,
            directional_budget_upper=p_upper, admissible_parameter_upper=t_cap,
            example_parameter=t_example).items()},
        example_parameter_dyadic_exponent=exponent,
        verdict='conditional-positive-same-action-budget-interval',
        conditions=list(dict.fromkeys(common['conditions']+primal['conditions']))+[
            'The same actual complete H-T identity and global prime/diagonal slack give H-T >= I/632.',
            'The selected Schur Pi is the same common actual form and inherits Pi >= I/632.',
            'The archived energies enclose the SAME unscaled full u=M w0, with constant Wfar=Dfar I.',
        ],
        scope='Conditional on the stated common actual forms/domains and sound archives, '
              'u u* <= p_upper Pi and every rational 0<t<=t_cap satisfies '
              'p_upper<=h(1-t)/t. This is an admissible directional shift budget only; '
              'no ell sign, sufficient update, actual retained/cofinal positivity, Robin-tail, '
              'RH, generic novelty or Lean claim.',
        old_scientific_productions_replayed=0, new_source_actions=0,
        new_inverse_or_solve_calls=0, parameter_scan=False,
        elapsed_seconds=time.monotonic()-started)
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('An input changed during consumption: '+name)
    with out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(result['verdict'], 'dyadic example 2^-'+str(exponent), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output; accepted packets are never overwritten.')
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
