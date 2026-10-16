"""Use the certified small reference inverse to contract shared error Grams.

Existing Neumann/inertia and Loewner comparison turn the cached exact
reference into trace bounds. Only thin products run: no inverse, solve,
scientific supplier, determinant or eigensolver is regenerated.
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
from full_gamma_directional_metric import balls, midpoint
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in
from structured_inverse_residual_interval import envelope


def consume(args):
    started = time.monotonic()
    certificate, certificate_hash = read_json(args.certificate)
    shared, shared_hash = read_json(args.shared)
    core, core_hash = read_json(args.core)
    cross, cross_hash = read_json(args.cross)
    if (certificate['schema'], shared['schema'], core['schema'], cross['schema'],
            certificate['verdict'], certificate['exact_reference_identities_checked']) != (
            'ACTUAL_C27_FULL_GAMMA_DIRECTIONAL_METRIC_V1', 'ACTUAL_C27_FULL_GAMMA_CORRELATED_SIGN_V1',
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1', 'ACTUAL_C27_FULL_GAMMA_DEFECT_CROSS_V1',
            'directional-center-enclosed', True):
        raise ValueError('Require the accepted same-core exact-reference certificate and shared errors.')
    if (certificate['packet_sha256']['shared'], certificate['packet_sha256']['core'],
            certificate['packet_sha256']['cross'], shared['packet_sha256']['core'],
            shared['packet_sha256']['cross'], certificate['inherited_packet_sha256']) != (
            shared_hash, core_hash, cross_hash, core_hash, cross_hash, shared['packet_sha256']):
        raise ValueError('Different exact reference, old core, center or shared Gram provenance.')
    expected = shared['input_sha256']
    if any(packet['input_sha256'] != expected for packet in (certificate, core, cross)):
        raise ValueError('Different scientific inputs.')
    for name, digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Scientific input hash mismatch: '+name)
    modes = shared['remaining_modes']
    if len(modes) != 18 or modes[0] != 0 or any(
            packet['remaining_modes'] != modes for packet in (certificate, core, cross)):
        raise ValueError('Require the same full old18-coordinate lift.')
    entries, bits = core['reduced_schur_intervals'], core['grid_bits']
    sigma = symmetric(matrix_in(entries, bits))
    n = len(modes)-1
    sigma0 = fmpq_mat([[midpoint(entries[min(i, j)][max(i, j)], bits)
                        for j in range(n+1)] for i in range(n+1)])
    a, a0 = sigma[0, 0], sigma0[0, 0]
    if not a < 0 or not a0 < 0:
        raise ArithmeticError('Require the accepted negative first pivot.')
    reference = fmpq_mat([[sigma0[i+1, j+1]-sigma0[i+1, 0]*sigma0[0, j+1]/a0
                          for j in range(n)] for i in range(n)])
    R = fmpq_mat([[fmpq(value) for value in row] for row in certificate['reference_inverse']])
    identity = fmpq_mat([[int(i == j) for j in range(n)] for i in range(n)])
    if R*reference != identity:
        raise ArithmeticError('The cached exact inverse does not match this reference matrix.')
    eta = interval_ball(certificate['scalar_bounds']['neumann_eta_upper'], certificate['grid_bits']).upper()
    if not 0 <= eta < 1:
        raise ArithmeticError('Require the certified outward Neumann upper strictly below one.')
    T = arb_mat([[-sigma[0, j+1]/a if i == 0 else arb(int(i == j+1))
                  for j in range(n)] for i in range(n+1)])
    Rballs = balls(R)
    traces = {}
    for name, key in (('finite', 'finite_error_gram_on_lift'), ('far', 'far_error_gram_on_lift')):
        C = symmetric(matrix_in(shared[key], shared['grid_bits']))
        positive_gram = symmetric(T.transpose()*C*T)
        RC = Rballs*positive_gram
        reference_trace = sum((RC[i, i] for i in range(n)), arb(0))
        if not reference_trace.is_finite() or not reference_trace.upper() >= 0:
            raise ArithmeticError('A nonnegative reference trace upper was not enclosed.')
        lower = (max(arb(0), reference_trace.lower())/(1+eta)).lower()
        upper = (reference_trace.upper()/(1-eta)).upper()
        traces[name] = dict(reference=reference_trace, lower=lower, upper=upper)
    center_energy = interval_ball(certificate['center_energy_interval'], certificate['grid_bits'])
    if not center_energy.lower() > 0:
        raise ArithmeticError('Require the accepted positive center-energy enclosure.')
    nc_lo, nc_up = center_energy.lower().sqrt().lower(), center_energy.upper().sqrt().upper()
    nf, nt = traces['finite']['upper'].sqrt().upper(), traces['far']['upper'].sqrt().upper()
    positive_lo = max(arb(0), (nc_lo-nf-nt).lower())**2
    positive_up = ((nc_up+nf+nt)**2).upper()
    b0 = interval_ball(cross['common_row_on_lift'][0][0], cross['grid_bits'])
    negative = b0**2/a
    qlo, qup = (negative+positive_lo).lower(), (negative+positive_up).upper()
    D = interval_ball(cross['common_denominator_upper'], cross['grid_bits'])
    slack = (-D.upper()-qup).lower()
    verdict = ('positive-good-metric-budget' if slack > 0 else
               'Dstar-target-excluded-even-at-zero-good-metric' if qlo+D.lower() > 0 else
               'inconclusive-reference-trace-sign')
    scalar = dict(neumann_eta_upper=eta, center_energy=center_energy,
                  finite_reference_trace=traces['finite']['reference'],
                  far_reference_trace=traces['far']['reference'],
                  finite_trace_lower=traces['finite']['lower'], finite_trace_upper=traces['finite']['upper'],
                  far_trace_lower=traces['far']['lower'], far_trace_upper=traces['far']['upper'],
                  positive_energy_lower=positive_lo, positive_energy_upper=positive_up,
                  negative_energy=negative, signed_inverse_metric_lower=qlo,
                  signed_inverse_metric_upper=qup, common_denominator=D,
                  Dstar_target_scalar_lower=(qlo+D.lower()).lower(),
                  Dstar_target_scalar_upper=(qup+D.upper()).upper(),
                  sufficient_good_metric_budget_lower=slack)
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_REFERENCE_TRACE_V1',
                  packet_sha256=dict(certificate=certificate_hash, shared=shared_hash,
                                     core=core_hash, cross=cross_hash),
                  input_sha256=expected, remaining_modes=modes,
                  grid_bits=args.output_bits, working_bits=args.precision,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scalar_bounds={name: envelope(value, args.output_bits) for name, value in scalar.items()},
                  signed_inverse_metric_interval=[envelope(qlo, args.output_bits)[0],
                                                   envelope(qup, args.output_bits)[1]],
                  verdict=verdict, elapsed_seconds=time.monotonic()-started,
                  conditions=certificate['conditions']+[
                      'Actual Pplus is positive and the reference path is nonsingular for eta<1; hence P0 is positive.',
                      'R/(1+eta) <= Pplus^-1 <= R/(1-eta) applies to the same actual PSD transported error Grams.'],
                  scope='Known inverse/Loewner trace comparison on the saved actual lift; no inverse or solve. '
                    'Positive-metric bounds and the full-row b0 interval are conservative simultaneous bounds. '
                    'They do not assert simultaneous attainment or actual extremal errors. Good-block Gamma, '
                    'full retained/cofinal positivity, Lean verification and RH remain unsupplied.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    for name, value in scalar.items():
        print(name, value.str(14), flush=True)
    print(verdict, 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(certificate='full-gamma-directional-metric.json.zlib',
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
