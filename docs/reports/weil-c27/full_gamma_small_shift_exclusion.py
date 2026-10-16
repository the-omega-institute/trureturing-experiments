"""Extend the cached fixed-trial lower over the saved surplus interval.

Only the new mixed scalar and its shifted far penalty are evaluated.
The old trial, lower, square-completion constants and eligibility budget
are read verbatim; no inverse, action, Gram or old producer is rerun.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time
import zlib

from flint import arb, ctx, fmpq

from diagonal_floor_witness import interval_ball
from full_gamma_surplus_shift_gate import minimum, rational
from structured_inverse_residual_interval import envelope


HERE = Path(__file__).resolve().parent
PINS = {
    'full-gamma-baseline-gap.json.zlib':
        'f7182de82c7f015252fc14cffc2b0e947ed6ff557f1f58de464dbb9ac88699bd',
    'full-gamma-actual-joint-lower.json.zlib':
        'f294428dc9685e04fca777a5dc71ac0374fc17a4c69246422e8b4325351da5a9',
    'full-gamma-lift-projection.json.zlib':
        'a41813b27cd4567af326d42ec77e3b8c8942474b80bdf5833f7bc05aff52f873',
    'full-gamma-projection-upper.json.zlib':
        'a0f6f2f222d35fa6a1c94c68d1573d9b7695eea03308974ed4481b03940993e3',
    'full-gamma-actual-channel-projection.json.zlib':
        '63e3c2e291940ca3549fa4a1303b55b4d92c08d594a50daf7ca2132087ad2745',
    'joint-gamma-6765-rho32-v2-schur.json':
        '13823fc56c3c06d79e5db97132afbdbd3fef79929ce59103436fc8dfae9166fd',
    'full-gamma-surplus-budget.json':
        '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f',
}
MODES = [0, 1, 2, 3, 4, 5, 6, 9, 10, 12, 14, 15, 18, 19, 24, 29, 33, 39]


def consume(directory, out):
    started = time.monotonic()
    packets = []
    for name, digest in PINS.items():
        raw = (directory/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted input: '+name)
        packets.append(json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw))
    gap, joint, lift, projection, actual, core, budget = packets
    if [p['schema'] for p in packets] != [
            'ACTUAL_C27_FULL_GAMMA_BASELINE_GAP_V1',
            'ACTUAL_C27_FULL_GAMMA_ACTUAL_JOINT_LOWER_V1',
            'ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
            'ACTUAL_C27_COMMON_PROJECTION_UPPER_V1',
            'ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1',
            'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1']:
        raise ValueError('Different accepted packet roles.')
    if (any(p['grid_bits'] != 640 or p['remaining_modes'] != MODES
            or p['input_sha256'] != gap['input_sha256'] for p in packets[:-1])
            or gap['verdict'] != 'same-Wbar-stronger-base-fixed-dual-update-excluded'
            or budget['verdict'] != 'conditional-positive-same-action-budget-interval'):
        raise ValueError('Different common source, lift, grids or accepted lower/budget.')
    for key, filename in dict(joint='full-gamma-actual-joint-lower.json.zlib',
                              projection='full-gamma-lift-projection.json.zlib',
                              actual='full-gamma-actual-channel-projection.json.zlib',
                              core='joint-gamma-6765-rho32-v2-schur.json').items():
        if gap['packet_sha256'][key] != PINS[filename]:
            raise ValueError('Different baseline dependency: '+key)
    if (projection['packet_sha256']['gap'] != PINS['full-gamma-baseline-gap.json.zlib']
            or projection['packet_sha256']['projection'] != PINS['full-gamma-lift-projection.json.zlib']
            or projection['packet_sha256']['core'] != PINS['joint-gamma-6765-rho32-v2-schur.json']
            or budget['packet_sha256']['full-gamma-projection-upper.json.zlib'] !=
            PINS['full-gamma-projection-upper.json.zlib']
            or budget['source_integer_sha256'] != projection['exact_source_integer_sha256'][0]
            or budget['source_bits'] != projection['exact_source_bits'][0]
            or budget['source_bits'] != 384 or budget['source_row'] != 0):
        raise ValueError('Different actual original unscaled source or budget correspondence.')
    for packet, field in ((lift, 'lift_gram'), (lift, 'far_gram_upper'),
                          (core, 'reduced_schur_intervals')):
        if len(packet[field]) != 18 or any(len(row) != 18 for row in packet[field]):
            raise ValueError('Incomplete same-lift matrix: '+field)
    b_rows = actual['actual_finite_row_on_lift']
    r_rows = projection['matrices']['source_on_lift']
    saved_x = joint['new_reference_profile']
    if (len(b_rows) != 1 or len(b_rows[0]) != 18 or len(r_rows[0]) != 18
            or len(saved_x) != 17 or any(len(row) != 1 for row in saved_x)):
        raise ValueError('Incomplete finite row, shift row or saved exact trial.')

    bounds = gap['scalar_bounds']
    a = interval_ball(bounds['shifted_negative_pivot'], 640)
    delta = interval_ball(bounds['lifted_gap_coefficient'], 640)
    a_low = fmpq(int(bounds['shifted_far_quadratic_coefficient_lower'][0]), 1 << 640)
    p0 = fmpq(int(bounds['shifted_completed_square_penalty_upper'][1]), 1 << 640)
    l0 = Fraction(int(bounds['stronger_baseline_core_joint_lower'][0]), 1 << 640)
    h_low = Fraction(*map(int, budget['exact_rationals']['h_lower']))
    tcap = Fraction(*map(int, budget['exact_rationals']['admissible_parameter_upper']))
    if not (a < 0 and delta > 0 and a_low > 0 and p0 >= 0
            and l0 > 0 and h_low > 0 and 0 < tcap < 1):
        raise ValueError('Require accepted strict pivot, square-completion and budget signs.')
    x = [arb(fmpq(row[0])) for row in saved_x]
    r = [interval_ball(pair, 640) for pair in r_rows[0]]
    b0 = interval_ball(b_rows[0][0], 640)
    # Only the first row of SigmaMax is needed. Xi is not substituted.
    c = [interval_ball(core['reduced_schur_intervals'][0][j], 640)
         +delta*interval_ball(lift['lift_gram'][0][j], 640) for j in range(1, 18)]
    gradient0 = (b0-sum((cj*xj for cj, xj in zip(c, x)), arb(0)))/a
    mixed = r[0]*gradient0+sum((rj*xj for rj, xj in zip(r[1:], x)), arb(0))
    curvature = r[0]**2/a
    g00 = interval_ball(lift['far_gram_upper'][0][0], 640)
    if not g00 >= 0:
        raise ValueError('Require the paid PSD far Gram diagonal enclosure.')
    kappa = arb(fmpq(33, 32))
    p1 = (kappa**2*(r[0]/a)**2*g00.upper()/arb(a_low)).upper()
    cross_penalty = (arb(p0)*p1).sqrt().upper()
    linear = arb(fmpq(h_low.numerator, h_low.denominator))+2*mixed-2*cross_penalty
    quadratic = curvature-p1
    scalar = dict(new_gradient_first_coordinate=gradient0, new_mixed_scalar=mixed,
                  new_signed_curvature=curvature, shifted_far_penalty_coefficient_upper=p1,
                  far_penalty_cross_sqrt_upper=cross_penalty,
                  lower_polynomial_linear=linear, lower_polynomial_quadratic=quadratic)
    exported = {key: envelope(value, 640) for key, value in scalar.items()}
    if any(value is None for value in exported.values()):
        raise ArithmeticError('Nonfinite new scalar enclosure; stop without a precision retry.')
    coefficients = [l0, Fraction(int(exported['lower_polynomial_linear'][0]), 1 << 640),
                    Fraction(int(exported['lower_polynomial_quadratic'][0]), 1 << 640)]
    unit_coefficients = [coefficients[0], coefficients[1]*tcap, coefficients[2]*tcap*tcap]
    unit_argmin, lower_minimum = minimum(unit_coefficients)
    verdict = ('conditional-small-surplus-interval-update-excluded' if lower_minimum > 0
               else 'inconclusive-small-surplus-interval-lower')
    result = dict(schema='ACTUAL_C27_SMALL_SHIFT_EXCLUSION_V1', packet_sha256=PINS,
                  input_sha256=gap['input_sha256'], remaining_modes=MODES,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  grid_bits=640, working_bits=1024, scalar_bounds=exported,
                  inherited_scalar_fields=dict(lower='gap.stronger_baseline_core_joint_lower',
                      far_penalty='gap.shifted_completed_square_penalty_upper',
                      far_quadratic='gap.shifted_far_quadratic_coefficient_lower',
                      h_lower='budget.h_lower', parameter_endpoint='budget.admissible_parameter_upper'),
                  coefficient_order=['constant', 'linear', 'quadratic'],
                  exact_lower_coefficients=list(map(rational, coefficients)),
                  exact_unit_interval_coefficients=list(map(rational, unit_coefficients)),
                  minimum_of_lower_polynomial=rational(lower_minimum),
                  minimizing_parameter=rational(unit_argmin*tcap),
                  parameter_interval=[rational(Fraction(0)), rational(tcap)],
                  verdict=verdict, exact_trial_source='joint.new_reference_profile',
                  old_scientific_productions_replayed=0, new_inverse_or_solve_calls=0,
                  parameter_scan=False, elapsed_seconds=time.monotonic()-started,
                  conditions=gap['conditions']+budget['conditions']+[
                      'The same real signed compression SigmaMax has exactly one negative eigenvalue; its positive complement, saved exact17trial, finite bJ, shift rV and complete paid PSD far Gram G share the same lift.',
                      'The full unscaled u and its far action stay fixed; K_tau=K+tau*h and B_tau=B+tau*w0*R_res on the same common domain.',
                      'The inherited completed-square constants and cached lower are sound, outward enclosures of this same fixed-trial argument; the metric triangle bounds the shifted far gradient.',
                      'SigmaMax upper-bounds compression of the actual same-Wbar stronger base; its inverse quadratic is used only for compression exclusion, not substituted into an arbitrary actual score.'],
                  scope='Conditional numerical lower for K_tau+qMax(B_tau V) on the saved sufficient '
                        'budget interval, including its zero endpoint. Strict positivity excludes '
                        'Lbar+B_tau*B_tau/D for every positive D>=K_tau in this family. '
                        'This is not actual Weil negativity or a maximal eligibility bound. '
                        'Larger admissible shifts, other sources/bases/majorants, retained/cofinal '
                        'positivity, Robin-tail, RH, generic novelty and Lean remain outside scope.')
    for name, digest in PINS.items():
        if hashlib.sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise ValueError('An input changed during consumption: '+name)
    with out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    for key, value in scalar.items():
        print(key, value.str(18), flush=True)
    print('lower-polynomial minimum', float(lower_minimum), verdict, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.json':
        parser.error('Use a fresh .json output; accepted packets are never overwritten.')
    ctx.prec = 1024
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
