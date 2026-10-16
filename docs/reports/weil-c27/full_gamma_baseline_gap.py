"""Test whether the paid far-enclosure gap can rescue the fixed dual.

Consume the accepted finite charge, saved exact trial, finite row and
lift/far Grams. Apply the existing variational inverse bound to the
changed eighteen-coordinate compression. No inverse, solve, old signed
energy, channel/action production or trace supplier is rerun.
"""
import argparse
import hashlib
from pathlib import Path
import time

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

from arb_inverse_residual_obstruction import HERE
from diagonal_floor_witness import interval_ball
from full_gamma_correlated_sign import symmetric
from full_gamma_directional_metric import balls
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, write_packet
from structured_inverse_residual_interval import envelope


def consume(args):
    started = time.monotonic()
    names = ('old', 'old_pivots', 'base', 'projection', 'actual', 'joint',
             'certificate', 'common', 'core')
    packets, hashes = {}, {}
    for name in names:
        packets[name], hashes[name] = read_json(getattr(args, name))
    old, pivots, base, projection, actual, joint, certificate, common, core = (
        packets[name] for name in names)
    schemas = dict(old='ACTUAL_C27_STRUCTURED_TARGET_V1',
                   base='ACTUAL_C27_STRUCTURED_TARGET_V1',
                   projection='ACTUAL_C27_FULL_GAMMA_LIFT_PROJECTION_V1',
                   actual='ACTUAL_C27_FULL_GAMMA_ACTUAL_CHANNEL_PROJECTION_V1',
                   joint='ACTUAL_C27_FULL_GAMMA_ACTUAL_JOINT_LOWER_V1',
                   certificate='ACTUAL_C27_FULL_GAMMA_DIRECTIONAL_METRIC_V1',
                   common='ACTUAL_C27_FULL_GAMMA_COMMON_COLUMN_V1',
                   core='ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1')
    if any(packets[name]['schema'] != value for name, value in schemas.items()):
        raise ValueError('Require accepted packets for the same OLD base and actual fixed dual.')
    if any(old['metadata'][name] != pivots[name] for name in (
            'far_matrix_error', 'complete_far_operator_allowance', 'trial_operator_norm_sq_upper')):
        raise ValueError('Original far allowance disagrees with its explicitly gridded packet.')
    if (pivots['grid_bits'], base['metadata']['original_metadata_grid_bits'],
            projection['original_metadata_grid_bits'], base['metadata']['common_young_parameter'],
            base['metadata']['original_target_sha256'], projection['original_target_sha256'],
            projection['base_target_sha256'], projection['base_core_sha256'],
            core['target_packet_sha256']) != (
            512, 512, 512, '1/32', hashes['old'], hashes['old'], hashes['base'],
            hashes['core'], hashes['base']):
        raise ValueError('Different base, core, original allowance grid or Young parameter.')
    if any(joint['packet_sha256'][name] != hashes[name] for name in (
            'actual', 'certificate', 'projection', 'common', 'core')):
        raise ValueError('Different accepted actual-source joint lower.')
    if (actual['packet_sha256']['projection'], actual['packet_sha256']['common'],
            certificate['packet_sha256']['core'], certificate['packet_sha256']['shared'],
            actual['packet_sha256']['shared'], joint['verdict'], certificate['verdict'],
            certificate['exact_reference_identities_checked']) != (
            hashes['projection'], hashes['common'], hashes['core'],
            joint['packet_sha256']['shared'], joint['packet_sha256']['shared'],
            'fixed-old-base-and-dual-rank-one-update-excluded', 'directional-center-enclosed', True):
        raise ValueError('Different reference or finite common source, or missing accepted sign evidence.')
    expected, modes = actual['input_sha256'], actual['remaining_modes']
    if len(modes) != 18 or modes[0] != 0 or any(
            packets[name]['input_sha256'] != expected or packets[name]['remaining_modes'] != modes
            for name in ('projection', 'joint', 'certificate', 'core')):
        raise ValueError('Different actual scientific inputs or exact18lift order.')
    if any(packets[name]['metadata']['input_sha256'] != expected for name in ('old', 'base', 'common')):
        raise ValueError('Different original target or common fixed dual.')
    for name, digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Scientific input hash mismatch: '+name)
    epsilon = fmpq(projection['exact_far_allowance'])
    original_epsilon = fmpq(int(old['metadata']['complete_far_operator_allowance'][1]), 1 << 512)
    eps_bits = base['metadata']['reused_far_allowance_grid_bits']
    base_epsilon = fmpq(int(base['metadata']['reused_far_allowance'][1]), 1 << eps_bits)
    if epsilon != base_epsilon or epsilon < original_epsilon or not epsilon > 0:
        raise ValueError('Require the same complete upward far allowance, not a midpoint or smaller grid.')

    n = len(modes)
    for key in ('lift_gram', 'far_gram_upper'):
        if len(projection[key]) != n or any(len(row) != n for row in projection[key]):
            raise ValueError('Incomplete same-lift Gram: '+key)
    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], core['grid_bits']))
    Hlift = symmetric(matrix_in(projection['lift_gram'], projection['grid_bits']))
    G = symmetric(matrix_in(projection['far_gram_upper'], projection['grid_bits']))
    kappa = arb(fmpq(33, 32))
    delta = 2*kappa*arb(epsilon)
    shifted = symmetric(sigma+delta*Hlift)
    a = shifted[0, 0]
    if (sigma.nrows(), sigma.ncols()) != (n, n) or not a < 0:
        raise ArithmeticError('Require the accepted old inertia and a strictly negative shifted pivot.')
    # The actual PSD addition preserves at least17positive eigenvalues.
    # The negative first quadratic value leaves exactly one negative;
    # hence this NEW complement is positive, without another inertia test.
    T = arb_mat([[-shifted[0, j+1]/a if i == 0 else arb(int(i == j+1))
                  for j in range(n-1)] for i in range(n)])
    P = symmetric(arb_mat([[shifted[i+1, j+1]-shifted[i+1, 0]*shifted[0, j+1]/a
                            for j in range(n-1)] for i in range(n-1)]))
    finite = matrix_in(actual['actual_finite_row_on_lift'], actual['grid_bits'])
    x = fmpq_mat([[fmpq(value) for value in row] for row in joint['new_reference_profile']])
    if (finite.nrows(), finite.ncols(), x.nrows(), x.ncols()) != (1, n, n-1, 1):
        raise ValueError('Require the actual finite row and the complete saved exact17trial.')
    xb = balls(x)
    h = finite*T
    trial_positive = 2*(h*xb)[0, 0]-(xb.transpose()*P*xb)[0, 0]
    qtrial = (finite[0, 0]**2/a+trial_positive).lower()
    gradient = T*xb
    gradient[0, 0] += finite[0, 0]/a
    gradient_energy = (gradient.transpose()*G*gradient)[0, 0]
    Bsq = (kappa**2*gradient_energy.upper()).upper()
    A = (kappa-kappa**2*G[0, 0].upper()/(-a).lower()).lower()
    if not A > 0 or not Bsq >= 0:
        raise ArithmeticError('The changed-source completed square needs A>0 and a nonnegative B-square ceiling.')
    penalty = (Bsq/A).upper()
    charge = interval_ball(joint['scalar_bounds']['actual_common_charge_floor'], joint['grid_bits']).lower()
    lower = (charge+qtrial-penalty).lower()
    scalar = dict(full_far_allowance=arb(epsilon), lifted_gap_coefficient=delta,
                  shifted_negative_pivot=a, shifted_positive_trial_value=trial_positive,
                  shifted_finite_signed_trial_lower=qtrial,
                  shifted_far_trial_gradient_energy=gradient_energy,
                  shifted_far_quadratic_coefficient_lower=A,
                  shifted_far_linear_coefficient_sq_upper=Bsq,
                  shifted_completed_square_penalty_upper=penalty,
                  inherited_common_charge_floor=charge, stronger_baseline_core_joint_lower=lower)
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_BASELINE_GAP_V1', packet_sha256=hashes,
                  input_sha256=expected, remaining_modes=modes, grid_bits=args.output_bits,
                  working_bits=args.precision, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scalar_bounds={name: envelope(value, args.output_bits) for name, value in scalar.items()},
                  verdict=('same-Wbar-stronger-base-fixed-dual-update-excluded' if lower > 0
                           else 'inconclusive-fixed-dual-baseline-gap'),
                  exact_trial_source='joint.new_reference_profile',
                  elapsed_seconds=time.monotonic()-started,
                  conditions=joint['conditions']+[
                      'The inherited two-sided far approximation has norm error at most the SAME upward epsilon, so 0<=Fup-Factual<=2*epsilon*I.',
                      'For these SAME Q,R,Wbar,Lstar, Delta=kappa*(Fhat-Factual)+(eps0+rho*eps1)*I with ||Fhat-Factual||<=eps0<=eps1<=epsilon, hence 0<=Delta<=2*kappa*epsilon*I. The actual original and transported paddings are not identified.',
                      'The accepted actual common charge floor bounds mJ-h-Ephi for the SAME fixed common source.',
                      'Exact old Sigma has inertia17positive/1negative. PSD addition of delta*Hlift preserves at least17positive eigenvalues, and the strictly negative shifted first pivot leaves exactly one negative.',
                      'The archived lift Gram is V*V, and the full signed far Gram plus its entire allowance is a PSD majorant on that SAME lift.',
                      'The saved exact17profile is only a trial for the CHANGED positive complement; its variational lower is used without asserting it is a new inverse action.',
                      'Actual t=kappa*ufar_adjoint*Wfar*Rfar*V and s=||ufar||W^2 share the SAME source; the displayed new gradient and G control their linear term.',
                      'For every positive D>=Kcal, strict positivity of D+b*SigmaMax^-1*b_adjoint contradicts PSD by the existing one-negative-core rank-one criterion.'],
                  scope='Conditional numerical exclusion of Lbar+calB*calB/D for Lbar=Q-R*WbarR '
                    'with exactly the old19channels/rho1/32 and saved fixeddyadicw, everypositiveD>=Kcal. '
                    'Other majorants, new or combineddual directions and actualS_N positivity remain outside '
                    'this exclusion. No inverse/solve/supplier rerun, generic novelty, Lean, Robin-tail or RH result.')
    write_packet(args.out, result)
    for name, value in scalar.items():
        print(name, value.str(18), flush=True)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(old='structured-joint-2048-target.json.zlib',
                             old_pivots='structured-joint-2048-pivots.json',
                             base='joint-gamma-6765-rho32-v2-target.json.zlib',
                             projection='full-gamma-lift-projection.json.zlib',
                             actual='full-gamma-actual-channel-projection.json.zlib',
                             joint='full-gamma-actual-joint-lower.json.zlib',
                             certificate='full-gamma-directional-metric.json.zlib',
                             common='full-gamma-common-column.json.zlib',
                             core='joint-gamma-6765-rho32-v2-schur.json').items():
        parser.add_argument('--'+name.replace('_', '-'), type=Path, default=HERE/default)
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
