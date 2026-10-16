"""Bound the fixed common-source joint scalar below using actual channels.

Reuse the exact17reference inverse for the newly supplied finite row.
Combine far energy with its same-source cross term before taking a
bound. No inverse, solve, old scalar consumer or scientific producer
runs. A positive result excludes only the specified old-base update.
"""
import argparse
import hashlib
from pathlib import Path
import time

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

from arb_inverse_residual_obstruction import HERE
from diagonal_floor_witness import interval_ball
from full_gamma_correlated_sign import symmetric
from full_gamma_directional_metric import balls, midpoint, rationals
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, write_packet
from structured_inverse_residual_interval import envelope


def consume(args):
    started = time.monotonic()
    packets, hashes = {}, {}
    for name in ('actual', 'certificate', 'trace', 'shared', 'core', 'projection', 'common', 'defect'):
        packets[name], hashes[name] = read_json(getattr(args, name))
    actual, certificate, trace, shared, core, projection, common, defect = (
        packets[name] for name in ('actual', 'certificate', 'trace', 'shared', 'core', 'projection', 'common', 'defect'))
    schemas = dict(actual='ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
                   certificate='ACTUAL_C27_FULL_GAMMA_DIRECTIONAL_METRIC_V1',
                   trace='ACTUAL_C27_FULL_GAMMA_REFERENCE_TRACE_V1',
                   shared='ACTUAL_C27_FULL_GAMMA_CORRELATED_SIGN_V1',
                   core='ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1',
                   projection='ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
                   common='ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1',
                   defect='ACTUAL_C27_FINITE_DUAL_DEFECT_V1')
    if any(packets[name]['schema'] != value for name, value in schemas.items()):
        raise ValueError('Require accepted actual channels and the matching old-core reference/trace.')
    if any(actual['packet_sha256'][name] != hashes[name] for name in (
            'shared', 'projection', 'common', 'defect')):
        raise ValueError('Different actual finite common source.')
    if (certificate['packet_sha256']['shared'], certificate['packet_sha256']['core'],
            trace['packet_sha256']['certificate'], trace['packet_sha256']['shared'],
            trace['packet_sha256']['core'], projection['base_core_sha256'],
            common['finite_defect_sha256'], shared['packet_sha256']['core'],
            certificate['verdict'], certificate['exact_reference_identities_checked']) != (
            hashes['shared'], hashes['core'], hashes['certificate'], hashes['shared'],
            hashes['core'], hashes['core'], hashes['defect'], hashes['core'],
            'directional-center-enclosed', True):
        raise ValueError('Different old positive complement, exact reference or shared far trace.')
    expected, modes = actual['input_sha256'], actual['remaining_modes']
    if len(modes) != 18 or modes[0] != 0 or any(
            packets[name]['input_sha256'] != expected or packets[name]['remaining_modes'] != modes
            for name in ('certificate', 'trace', 'shared', 'core', 'projection')):
        raise ValueError('Different source or old18coordinate lift.')
    if common['metadata'] != defect['metadata'] or common['metadata']['input_sha256'] != expected:
        raise ValueError('Different actual middle energy and dual denominator.')
    for name, digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Scientific input hash mismatch: '+name)
    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], core['grid_bits']))
    a = sigma[0, 0]
    if not a < 0:
        raise ArithmeticError('Require the accepted negative first pivot.')
    n = len(modes)-1
    T = arb_mat([[-sigma[0, j+1]/a if i == 0 else arb(int(i == j+1))
                  for j in range(n)] for i in range(n+1)])
    P = symmetric(arb_mat([[sigma[i+1, j+1]-sigma[i+1, 0]*sigma[0, j+1]/a
                            for j in range(n)] for i in range(n)]))
    R = fmpq_mat([[fmpq(value) for value in row] for row in certificate['reference_inverse']])
    if (R.nrows(), R.ncols()) != (n, n):
        raise ValueError('Incomplete accepted exact reference inverse.')
    etaN = interval_ball(certificate['scalar_bounds']['neumann_eta_upper'], certificate['grid_bits']).upper()
    if not 0 <= etaN < 1:
        raise ArithmeticError('Require the accepted Neumann upper strictly below one.')
    finite = matrix_in(actual['actual_finite_row_on_lift'], actual['grid_bits'])
    if (finite.nrows(), finite.ncols()) != (1, n+1):
        raise ValueError('Incomplete actual finite common row.')
    h = finite*T
    h0 = fmpq_mat([[midpoint(envelope(h[0, j], args.output_bits), args.output_bits)
                    for j in range(n)]])
    x = R*h0.transpose()
    xb, Rb = balls(x), balls(R)
    residual = h.transpose()-P*xb
    residual_one = sum((abs(residual[i, 0]).upper() for i in range(n)), arb(0)).upper()
    Rres = Rb*residual
    residual_infinity = max(abs(Rres[i, 0]).upper() for i in range(n))
    Qc = (residual_one*residual_infinity/(1-etaN)).upper()
    variational = 2*(h*xb)[0, 0]-(xb.transpose()*P*xb)[0, 0]
    positive_lower = max(arb(0), variational.lower())
    positive_upper = (variational.upper()+Qc).upper()
    qfinite = finite[0, 0]**2/a+arb(positive_lower, 0)
    qlo = qfinite.lower()
    qup = (finite[0, 0]**2/a+positive_upper).upper()
    gradient = T*xb
    gradient[0, 0] += finite[0, 0]/a
    Ct = symmetric(matrix_in(shared['far_error_gram_on_lift'], shared['grid_bits']))
    trial_metric = (gradient.transpose()*Ct*gradient)[0, 0]
    tau = interval_ball(trace['scalar_bounds']['far_trace_upper'], trace['grid_bits']).upper()
    if not trial_metric.is_finite() or not trial_metric.upper() >= 0 or not tau >= 0:
        raise ArithmeticError('Require finite nonnegative actual far-seminorm ceilings.')
    gradient_remainder = (tau*Qc).sqrt().upper()
    Lfar = (trial_metric.upper().sqrt().upper()+gradient_remainder).upper()
    EMlo = interval_ball(common['bounds']['far_action_energy_upper'], common['grid_bits']).lower()
    if not EMlo > 0:
        raise ArithmeticError('Require a positive lower for the same fixed far-energy cap used in Ct.')
    G00 = interval_ball(projection['far_gram_upper'][0][0], projection['grid_bits']).upper()
    kappa = arb(fmpq(33, 32))
    A_lower = (kappa-kappa**2*G00/(-a).lower()).lower()
    if not A_lower > 0:
        raise ArithmeticError('This completed-square lower bound requires a positive far coefficient.')
    Bsq_upper = (Lfar**2/EMlo).upper()
    far_penalty = (Bsq_upper/A_lower).upper()
    middle = interval_ball(common['bounds']['middle_action_energy'], common['grid_bits'])
    dual = interval_ball(defect['denominator'], defect['grid_bits'])
    Ephi = interval_ball(actual['scalar_bounds']['actual_channel_energy'], actual['grid_bits'])
    charge_floor = (middle.lower()-dual.upper()-Ephi.upper()).lower()
    center_floor = (charge_floor+qlo).lower()
    joint_lower = (center_floor-far_penalty).lower()
    scalar = dict(neumann_eta_upper=etaN, new_positive_energy_lower=positive_lower,
                  new_positive_energy_upper=positive_upper, new_inverse_residual_charge_upper=Qc,
                  actual_finite_signed_energy_lower=qlo, actual_finite_signed_energy_upper=qup,
                  new_far_trial_gradient_energy=trial_metric,
                  new_far_gradient_remainder_upper=gradient_remainder,
                  new_far_gradient_norm_upper=Lfar, far_energy_cap_lower=EMlo,
                  far_quadratic_coefficient_lower=A_lower, far_linear_coefficient_sq_upper=Bsq_upper,
                  completed_square_far_penalty_upper=far_penalty,
                  actual_common_charge_floor=charge_floor,
                  actual_finite_joint_floor=center_floor, actual_joint_scalar_lower=joint_lower)
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_ACTUAL_JOINT_LOWER_V1',
                  packet_sha256=hashes, input_sha256=expected, remaining_modes=modes,
                  grid_bits=args.output_bits, working_bits=args.precision,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  new_reference_profile=rationals(x), new_reference_rhs=rationals(h0.transpose()),
                  scalar_bounds={name: envelope(value, args.output_bits) for name, value in scalar.items()},
                  verdict=('fixed-old-base-and-dual-rank-one-update-excluded' if joint_lower > 0
                           else 'inconclusive-actual-common-source-joint-lower'),
                  elapsed_seconds=time.monotonic()-started,
                  conditions=actual['conditions']+trace['conditions']+[
                      'The SAME full far Gram including its signed core and entire allowance controls the actual R_far V.',
                      'Ct=kappa^2*EM0*G for the SAME fixed cap EM0 from the accepted common packet; EM0 is at least the saved positive EMlo.',
                      'Actual t=kappa*u_far*W_far*R_far V and s=||u_far||_W^2 share the SAME source, so |t0|^2<=kappa^2*s*G00 and |t g|^2<=kappa^2*s*g*Gg.',
                      'The actual common charge is middle_action_energy-h+kappa*s-Ephi, and Gamma>=0.',
                      'Only the fixed OLD base Lstar and fixed saved dyadic dual w are tested; a positive result excludes PSD of Lstar+calB*calB/D for every positive D>=Kcal via the accepted one-negative-core criterion.'],
                  scope='A lower bound for the actual fixed-source joint scalar Kcal+Gamma+q '
                    'under inherited common-form/domain/inverse/comparison/sound-cache premises. '
                    'A positive verdict excludes only the specified fixed-old-base/fixed-dual rank-one '
                    'sufficient lower construction. It does not exclude the stronger Q-R*WbarR base, '
                    'other dual directions or general rank-one methods, and is not actual S_N '
                    'negativity, full retained/cofinal positivity, Robin-tail or RH. '
                    'Existing residual and weighted square-completion methods; no generic novelty or Lean claim.')
    write_packet(args.out, result)
    for name, value in scalar.items():
        print(name, value.str(18), flush=True)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(actual='full-gamma-actual-channel-projection.json.zlib',
                             certificate='full-gamma-directional-metric.json.zlib',
                             trace='full-gamma-reference-trace.json',
                             shared='full-gamma-correlated-sign.json',
                             core='joint-gamma-6765-rho32-v2-schur.json',
                             projection='full-gamma-lift-projection.json.zlib',
                             common='full-gamma-common-column.json.zlib',
                             defect='full-gamma-dual-32768-defect.json.zlib').items():
        parser.add_argument('--'+name, type=Path, default=HERE/default)
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib' or args.output_bits != 640 or args.precision < 832:
        parser.error('Use a fresh .zlib output and the inherited640-bit grid with192guardbits.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
