"""Use shared error correlations in the signed inverse quadratic form.

Reuse the certified center profile, residual energy, error traces and
full error Grams. Quadratic expansion and weighted Cauchy bound the
joint signed scalar without an inverse, solve or scientific supplier.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

from arb_inverse_residual_obstruction import HERE
from diagonal_floor_witness import interval_ball
from full_gamma_correlated_sign import symmetric
from full_gamma_directional_metric import balls
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in
from structured_inverse_residual_interval import envelope


def consume(args):
    started = time.monotonic()
    packets, hashes = {}, {}
    for name in ('certificate', 'trace', 'shared', 'core', 'cross'):
        packets[name], hashes[name] = read_json(getattr(args, name))
    certificate, trace, shared, core, cross = (packets[name] for name in (
        'certificate', 'trace', 'shared', 'core', 'cross'))
    schemas = dict(certificate='ACTUAL_C27_FULL_GAMMA_DIRECTIONAL_METRIC_V1',
                   trace='ACTUAL_C27_FULL_GAMMA_REFERENCE_TRACE_V1',
                   shared='ACTUAL_C27_FULL_GAMMA_CORRELATED_SIGN_V1',
                   core='ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1',
                   cross='ACTUAL_C27_FULL_GAMMA_DEFECT_CROSS_V1')
    if any(packets[name]['schema'] != value for name, value in schemas.items()):
        raise ValueError('Require accepted center, trace, shared Gram and old-core packets.')
    for name in ('certificate', 'shared', 'core', 'cross'):
        if trace['packet_sha256'][name] != hashes[name]:
            raise ValueError('Different reference-trace provenance: '+name)
    if (certificate['packet_sha256']['shared'], certificate['packet_sha256']['core'],
            certificate['packet_sha256']['cross'], shared['packet_sha256']['core'],
            shared['packet_sha256']['cross'], certificate['verdict']) != (
            hashes['shared'], hashes['core'], hashes['cross'], hashes['core'], hashes['cross'],
            'directional-center-enclosed'):
        raise ValueError('Different exact profile, same-object center or shared error Grams.')
    expected, modes = shared['input_sha256'], shared['remaining_modes']
    if len(modes) != 18 or modes[0] != 0 or any(
            packet['remaining_modes'] != modes or packet['input_sha256'] != expected
            for packet in (certificate, trace, core, cross)):
        raise ValueError('Different old-lift coordinates or scientific inputs.')
    if [row['mode'] for row in cross['rows']] != modes:
        raise ValueError('Different saved finite-center order.')
    for name, digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Scientific input hash mismatch: '+name)
    sigma = symmetric(matrix_in(core['reduced_schur_intervals'], core['grid_bits']))
    a = sigma[0, 0]
    if not a < 0:
        raise ArithmeticError('Require the accepted negative pivot.')
    n = len(modes)-1
    T = arb_mat([[-sigma[0, j+1]/a if i == 0 else arb(int(i == j+1))
                  for j in range(n)] for i in range(n+1)])
    x = balls(fmpq_mat([[fmpq(value) for value in row] for row in certificate['reference_profile']]))
    if x.nrows() != n or x.ncols() != 1:
        raise ValueError('Incomplete accepted exact reference profile.')
    bh0 = interval_ball(cross['rows'][0]['finite_center'], cross['grid_bits'])
    gradient_trial = T*x
    gradient_trial[0, 0] += bh0/a
    Qc = interval_ball(certificate['scalar_bounds']['inverse_residual_charge_upper'],
                       certificate['grid_bits']).upper()
    if not Qc >= 0:
        raise ArithmeticError('Require the accepted nonnegative residual-energy ceiling.')
    center_energy = interval_ball(certificate['center_energy_interval'], certificate['grid_bits'])
    center_q = bh0**2/a+center_energy
    scalar = dict(center_signed_energy=center_q, inverse_residual_charge_upper=Qc)
    linear, zero_radii, positive_radii = [], [], []
    for name, key in (('finite', 'finite_error_gram_on_lift'), ('far', 'far_error_gram_on_lift')):
        C = symmetric(matrix_in(shared[key], shared['grid_bits']))
        trial_energy = (gradient_trial.transpose()*C*gradient_trial)[0, 0]
        tau = interval_ball(trace['scalar_bounds'][name+'_trace_upper'], trace['grid_bits']).upper()
        if not trial_energy.is_finite() or not trial_energy.upper() >= 0 or not tau >= 0 or not C[0, 0].upper() >= 0:
            raise ArithmeticError('Require finite nonnegative seminorm and trace ceilings.')
        remainder = (tau*Qc).sqrt().upper()
        bound = (trial_energy.upper().sqrt().upper()+remainder).upper()
        linear.append(bound)
        zero_radii.append(C[0, 0].upper().sqrt().upper())
        positive_radii.append(tau.sqrt().upper())
        scalar[name+'_trial_gradient_energy'] = trial_energy
        scalar[name+'_gradient_residual_norm_upper'] = remainder
        scalar[name+'_linear_error_abs_upper'] = bound
    linear_upper = (2*sum(linear, arb(0))).upper()
    quadratic_negative = (sum(zero_radii, arb(0))**2/(-a).lower()).upper()
    quadratic_positive = (sum(positive_radii, arb(0))**2).upper()
    qlo = (center_q.lower()-linear_upper-quadratic_negative).lower()
    qup = (center_q.upper()+linear_upper+quadratic_positive).upper()
    D = interval_ball(cross['common_denominator_upper'], cross['grid_bits'])
    if not D.lower() > 0:
        raise ArithmeticError('Require the accepted positive common-denominator enclosure.')
    target_lo, target_up = (D.lower()+qlo).lower(), (D.upper()+qup).upper()
    slack = (-D.upper()-qup).lower()
    verdict = ('positive-good-metric-budget' if slack > 0 else
               'Dstar-target-excluded-even-at-zero-good-metric' if target_lo > 0 else
               'inconclusive-joint-signed-bound')
    scalar.update(linear_error_abs_upper=linear_upper,
                  quadratic_negative_abs_upper=quadratic_negative,
                  quadratic_positive_upper=quadratic_positive,
                  signed_inverse_metric_lower=qlo, signed_inverse_metric_upper=qup,
                  common_denominator=D, Dstar_target_scalar_lower=target_lo,
                  Dstar_target_scalar_upper=target_up,
                  sufficient_good_metric_budget_lower=slack)
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_JOINT_SIGN_V1', packet_sha256=hashes,
                  input_sha256=expected, remaining_modes=modes,
                  grid_bits=args.output_bits, working_bits=args.precision,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scalar_bounds={name: envelope(value, args.output_bits) for name, value in scalar.items()},
                  signed_inverse_metric_interval=[envelope(qlo, args.output_bits)[0],
                                                   envelope(qup, args.output_bits)[1]],
                  verdict=verdict, elapsed_seconds=time.monotonic()-started,
                  conditions=trace['conditions']+[
                      'The center residual is h*-Pplus*x for the SAME h and Pplus as the accepted center certificate.',
                      'Both full error Grams control components of the SAME actual common row, including negative/positive cross correlations.',
                      'The signed expansion uses uniform bounds, without asserting joint attainment of separate error extrema.'],
                  scope='Joint signed quadratic bounds for the particular Dstar-based sufficient lower update. '
                    'Strict positive Dstar+q lower excludes that update for every Gamma>=0; it does not '
                    'give actual Schur negativity, exclude every rank-one construction or refute RH. '
                    'A positive good-metric budget would still require Gamma to be supplied. '
                    'No generic novelty, Lean or full retained/cofinal positivity claim.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    for name, value in scalar.items():
        print(name, value.str(14), flush=True)
    print(verdict, 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(certificate='full-gamma-directional-metric.json.zlib',
                             trace='full-gamma-reference-trace.json',
                             shared='full-gamma-correlated-sign.json',
                             core='joint-gamma-6765-rho32-v2-schur.json',
                             cross='full-gamma-defect-cross.json').items():
        parser.add_argument('--'+name, type=Path, default=HERE/default)
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.output_bits != 640 or args.precision < 832:
        parser.error('Use a fresh output and the inherited640-bit grid with192 guard bits.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
