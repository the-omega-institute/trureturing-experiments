"""Extend the saved fixed-charge exclusion using the accepted charge floor.

Only exact rational scalar arithmetic is new. All common-source estimates,
reference errors, projection floors and matrix calculations are reused.
Run with Python 3.10+ and --out pointing to a fresh .json.zlib file.
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
    'full-gamma-conditioned-norm-budget.json.zlib': '5e0cbe77aac75340a29b40987bf820efa3272e6b421c40ba769cc38cdd9aef9a',
    'full-gamma-far-energy-budget.json.zlib': '11ccabd6110436b2b29f1a910c6c36a5fd8479af19d183cf614826ed9ad746c6',
    'full-gamma-actual-joint-lower.json.zlib': 'f294428dc9685e04fca777a5dc71ac0374fc17a4c69246422e8b4325351da5a9',
    'full-gamma-actual-channel-projection.json.zlib': '63e3c2e291940ca3549fa4a1303b55b4d92c08d594a50daf7ca2132087ad2745',
    'joint-gamma-6765-rho32-v2-schur.json': '13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd',
}
NORM, BUDGET, JOINT, ACTUAL, CORE = PINS
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'


def rational(packet, key):
    pair = packet['exact_rationals'][key]
    return Fraction(int(pair[0]), int(pair[1]))


def consume(out):
    started = time.perf_counter()
    packets = {}
    for name, digest in PINS.items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Saved input hash mismatch: '+name)
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    norm, budget, joint, actual, core = packets.values()
    schemas = ['ACTUAL_C27_CONDITIONED_SCALAR_NORM_BUDGET_V1',
               'ACTUAL_C27_FAR_INVERSE_ENERGY_BUDGET_V1',
               'ACTUAL_C27_FULL_GAMMA_ACTUAL_JOINT_LOWER_V1',
               'ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
               'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1']
    if ([p['schema'] for p in packets.values()] != schemas
            or any(p['grid_bits'] != 640 for p in packets.values())
            or len(core['remaining_modes']) != 18
            or any(p['remaining_modes'] != core['remaining_modes'] for p in packets.values())
            or any((p['source_integer_sha256'], p['source_bits'], p['source_row'])
                   != (SOURCE_ID, 384, 0) for p in (norm, budget))):
        raise ValueError('Different accepted comparison, grid, source or lift order.')
    if (any(norm['packet_sha256'][k] != PINS[k] for k in (BUDGET, ACTUAL, CORE))
            or any(budget['packet_sha256'][k] != PINS[k] for k in (ACTUAL, CORE))
            or joint['packet_sha256']['actual'] != PINS[ACTUAL]
            or joint['packet_sha256']['core'] != PINS[CORE]
            or joint['input_sha256'] != actual['input_sha256']
            or actual['input_sha256'] != core['input_sha256']
            or norm['verdict'] != 'fixed_optimistic_comparison_excluded_by_scalar_norm'
            or joint['verdict'] != 'fixed-old-base-and-dual-rank-one-update-excluded'):
        raise ValueError('Missing SAME-source accepted charge-floor/comparison chain.')
    old_upper = rational(norm, 'complete_joint_budget_upper')
    charge_upper = rational(budget, 'fixed_charge')
    charge_floor = Fraction(int(joint['scalar_bounds']['actual_common_charge_floor'][0]), 1 << 640)
    t_plus = rational(budget, 'comparison_t_plus')
    a_lower = rational(norm, 'actual_comparison_denominator_lower')
    if not (old_upper < 0 and 0 < charge_floor <= charge_upper and t_plus > 0 and a_lower > 0):
        raise ArithmeticError('Require the saved negative upper and positive comparison domains.')
    # zcrit(D)-QJ = zcrit(Dbar)-QJ + Dbar-D, with all other data fixed.
    # K >= charge_floor and D >= K, so this applies to every valid charge.
    charge_debit = charge_upper-charge_floor
    upper = old_upper+charge_debit
    exact = dict(inherited_fixed_charge_upper=charge_upper,
                 inherited_actual_common_charge_floor=charge_floor,
                 inherited_complete_budget_upper=old_upper,
                 charge_reduction_debit=charge_debit,
                 complete_admissible_charge_budget_upper=upper,
                 comparison_t_plus=t_plus,
                 comparison_denominator_lower=a_lower)
    if upper < 0:
        exact['endpoint_comparison_determinant_positive_lower'] = a_lower*(-upper)
    result = dict(schema='ACTUAL_C27_ADMISSIBLE_CHARGE_EXCLUSION_V1',
        packet_sha256=PINS, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source_integer_sha256=SOURCE_ID, source_bits=384, source_row=0,
        remaining_modes=core['remaining_modes'], inherited_grid_bits=640,
        exact_rationals={k: [str(v.numerator), str(v.denominator)] for k, v in exact.items()},
        verdict=('same_old_base_all_valid_charges_excluded'
                 if upper < 0 else 'inconclusive_admissible_charge_extension'),
        parameter_scope='same source and old lower target; 0<t<=t_plus and D>=K>=charge_floor',
        new_source_far_target_inverse_solve_queries=0,
        new_vector_matrix_reference_actions_or_contractions=0,
        old_scientific_producers_or_fixtures_replayed=0,
        source_direction_parameter_precision_cutoff_scan=False,
        comparison_coefficient_eligibility_claim=False,
        actual_full_border_sign_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started,
        conditions=[
            'All inherited common-form/domain, physical-majorant, source/lift, reference, complete-far and actual-denominator premises of the accepted norm comparison and charge floor remain required.',
            'The SAME actual charge obeys K=mJ-h+kappa*s-Ephi>=charge_floor, with s>=0. This is a lower bound, not the refined charge upper.',
            'At t_plus the complete saved upper changes by exactly Dbar-D; no paid estimate or nominal inverse is recomputed. The existing projection floor z>=QJ retains its lower-bound direction.',
            'If the new upper is negative, at D=charge_floor the endpoint border has positive first entry and determinant, hence is positive definite. The inherited inertia identity then gives a negative direction of the old lower target.',
            'Loewner order transports that negative direction to every 0<t<=t_plus and D>=K>=charge_floor on the SAME old lower target. The optimistic endpoint need not itself be eligible.',
            'No stronger true baseline, other source, larger shift, actual retained-form sign, cofinal positivity, Robin or RH result follows. This applies existing affine dependence and monotonicity, with no generic novelty or Lean claim.'])
    out.write_bytes(zlib.compress(json.dumps(result, sort_keys=True, separators=(',', ':')).encode(), 9))
    print(result['verdict'])
    print('bytes', out.stat().st_size, 'sha256', hashlib.sha256(out.read_bytes()).hexdigest())
    print('seconds', result['elapsed_seconds'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or not args.out.name.endswith('.json.zlib'):
        parser.error('Require a fresh .json.zlib output; no overwrite or retry.')
    consume(args.out)


if __name__ == '__main__':
    main()
