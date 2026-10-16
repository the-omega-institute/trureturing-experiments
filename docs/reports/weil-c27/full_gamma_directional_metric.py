"""Certify the center metric with an exact small reference and its residual.

Apply existing inverse-residual and Neumann bounds to the saved center
and old core. FLINT supplies the exact 17-dimensional reference inverse;
no old solve, scientific supplier, determinant or eigensolver is replayed.
"""
import argparse
import hashlib
from pathlib import Path
import time

from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

from arb_inverse_residual_obstruction import HERE
from diagonal_floor_witness import interval_ball
from full_gamma_correlated_sign import symmetric
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, write_packet
from structured_inverse_residual_interval import envelope


def midpoint(pair, bits):
    return fmpq(int(pair[0])+int(pair[1]), 1 << (bits+1))


def balls(matrix):
    return arb_mat([[arb(matrix[i, j]) for j in range(matrix.ncols())]
                    for i in range(matrix.nrows())])


def rationals(matrix):
    return [[str(matrix[i, j]) for j in range(matrix.ncols())]
            for i in range(matrix.nrows())]


def consume(args):
    started = time.monotonic()
    core, core_hash = read_json(args.core)
    cross, cross_hash = read_json(args.cross)
    rank, rank_hash = read_json(args.rank)
    shared, shared_hash = read_json(args.shared)
    if (core['schema'], cross['schema'], rank['schema'], shared['schema']) != (
            'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1', 'ACTUAL_C27_FULL_GAMMA_DEFECT_CROSS_V1',
            'ACTUAL_C27_SMALL_SCHUR_RANK_FLOOR_V1', 'ACTUAL_C27_FULL_GAMMA_CORRELATED_SIGN_V1'):
        raise ValueError('Require accepted old-core, cached center and shared-Gram packets.')
    parent = shared['packet_sha256']
    if (parent['core'], parent['cross'], parent['rank'], cross['base_core_sha256'],
            rank['schur_result_sha256'], rank['target_packet_sha256'],
            rank['pivot_result_sha256'], core['target_packet_sha256']) != (
            core_hash, cross_hash, rank_hash, core_hash, core_hash,
            cross['base_target_sha256'], core['pivot_result_sha256'], cross['base_target_sha256']):
        raise ValueError('Different old base, exact lift, center or inertia evidence.')
    paths = dict(core=args.core, cross=args.cross, rank=args.rank,
                 projection=HERE/'full-gamma-lift-projection.json.zlib',
                 common=HERE/'full-gamma-common-column.json.zlib',
                 boundary=HERE/'joint-gamma-6765-finite.json.zlib',
                 base=HERE/'joint-gamma-6765-rho32-v2-target.json.zlib',
                 defect=HERE/'full-gamma-dual-32768-defect.json.zlib',
                 direction=HERE/'joint-gamma-6765-direction.json')
    if set(paths) != set(parent):
        raise ValueError('Incomplete inherited common-object provenance.')
    for name, path in paths.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != parent[name]:
            raise ValueError('Inherited packet hash mismatch: '+name)
    expected = shared['input_sha256']
    if core['input_sha256'] != expected or cross['input_sha256'] != expected:
        raise ValueError('Different scientific inputs.')
    for name, digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Scientific input hash mismatch: '+name)
    modes = shared['remaining_modes']
    if any(packet['remaining_modes'] != modes for packet in (core, cross, rank)):
        raise ValueError('Different exact old-lift coordinate order.')
    if (len(modes), modes[0], rank['verdict'], rank['negative_count_lower'],
            rank['reused_positive_coordinate_count'], len(rank['leading_minors']), cross['rho']) != (
            18, 0, 'full-complement-inertia-counted', 1, 1263, 18, '1/32'):
        raise ValueError('Require the accepted old18-core with one negative direction.')
    for k, minor in enumerate(rank['leading_minors'], 1):
        if (minor['order'], minor['determinant_sign'], minor['pivot_sign']) != (
                k, -1, -1 if k == 1 else 1) or int(minor['determinant_interval'][1]) >= 0:
            raise ValueError('The accepted Jacobi evidence does not supply the positive complement.')
    if [row['mode'] for row in cross['rows']] != modes:
        raise ValueError('The saved finite centers do not match the lift order.')
    n = len(modes)-1
    entries, bits = core['reduced_schur_intervals'], core['grid_bits']
    sigma = symmetric(matrix_in(entries, bits))
    sigma0 = fmpq_mat([[midpoint(entries[min(i, j)][max(i, j)], bits)
                        for j in range(n+1)] for i in range(n+1)])
    a, a0 = sigma[0, 0], sigma0[0, 0]
    if not a < 0 or not a0 < 0:
        raise ArithmeticError('Require the accepted negative pivot and nonzero reference pivot.')
    T = arb_mat([[-sigma[0, j+1]/a if i == 0 else arb(int(i == j+1))
                  for j in range(n)] for i in range(n+1)])
    T0 = fmpq_mat([[-sigma0[0, j+1]/a0 if i == 0 else int(i == j+1)
                   for j in range(n)] for i in range(n+1)])
    positive = symmetric(arb_mat([[sigma[i+1, j+1]-sigma[i+1, 0]*sigma[0, j+1]/a
                                   for j in range(n)] for i in range(n)]))
    reference = fmpq_mat([[sigma0[i+1, j+1]-sigma0[i+1, 0]*sigma0[0, j+1]/a0
                          for j in range(n)] for i in range(n)])
    center_pairs = [row['finite_center'] for row in cross['rows']]
    center = arb_mat([[interval_ball(pair, cross['grid_bits']) for pair in center_pairs]])
    center0 = fmpq_mat([[midpoint(pair, cross['grid_bits']) for pair in center_pairs]])
    h, h0 = center*T, center0*T0
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_DIRECTIONAL_METRIC_V1',
                  packet_sha256=dict(core=core_hash, cross=cross_hash, rank=rank_hash, shared=shared_hash),
                  inherited_packet_sha256=parent, input_sha256=expected, remaining_modes=modes,
                  positive_dimension=n, working_bits=args.precision, grid_bits=args.output_bits,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  conditions=list(dict.fromkeys(shared['conditions']))+[
                      'The accepted actual Pplus is positive definite; a rational midpoint is only a reference.',
                      'The saved finite_center fields enclose the same exact bhat used by the shared-Gram consumer.',
                      'The Neumann residual bound applies only when the outward eta upper is strictly below one.'],
                  scope='A directional certificate for the center metric only, with measured input/residual '
                    'charges for this particular enclosure. No actual condition number, intrinsic uncertainty '
                    'lower bound, full-row error-trace sign, good-block metric, Lean or RH claim.')
    try:
        R = reference.inv()
        identity = fmpq_mat([[int(i == j) for j in range(n)] for i in range(n)])
        x = R*h0.transpose()
        if R*reference != identity or reference*x != h0.transpose():
            raise ArithmeticError('Exact rational inverse/profile identities failed.')
        k0 = (h0*x)[0, 0]
        result.update(reference_inverse=rationals(R), reference_profile=rationals(x),
                      reference_rhs=rationals(h0.transpose()), reference_energy=str(k0),
                      exact_reference_identities_checked=True)
        dp, dh = positive-balls(reference), h-balls(h0)
        E = arb_mat([[abs(dp[i, j]).upper() for j in range(n)] for i in range(n)])
        Rballs, xballs = balls(R), balls(x)
        Rabs = arb_mat([[abs(Rballs[i, j]).upper() for j in range(n)] for i in range(n)])
        K = Rabs*E
        eta = max(sum((K[i, j] for j in range(n)), arb(0)).upper() for i in range(n))
        scalar = dict(neumann_eta_upper=eta, reference_energy=arb(k0),
                      reference_profile_max_abs_upper=max(abs(xballs[i, 0]).upper() for i in range(n)),
                      reference_inverse_infinity_upper=max(sum((Rabs[i, j] for j in range(n)), arb(0)).upper()
                                                           for i in range(n)))
        if eta < 1:
            residual = dh.transpose()-dp*xballs
            one_norm = sum((abs(residual[i, 0]).upper() for i in range(n)), arb(0)).upper()
            Rr = Rballs*residual
            Rr_norm = max(abs(Rr[i, 0]).upper() for i in range(n))
            debit = (one_norm*Rr_norm/(1-eta)).upper()
            linear = 2*(dh*xballs)[0, 0]
            matrix = (xballs.transpose()*dp*xballs)[0, 0]
            base = arb(k0)+linear-matrix
            lower, upper = max(arb(0), base.lower()), (base.upper()+debit).upper()
            result['center_energy_interval'] = [envelope(lower, args.output_bits)[0],
                                                envelope(upper, args.output_bits)[1]]
            scalar.update(linear_perturbation=linear, matrix_perturbation=matrix,
                          variational_base=base, residual_one_norm_upper=one_norm,
                          preconditioned_residual_infinity_upper=Rr_norm,
                          inverse_residual_charge_upper=debit,
                          center_energy_lower=lower, center_energy_upper=upper)
            result['verdict'] = 'directional-center-enclosed'
        else:
            result['verdict'] = 'inconclusive-neumann-regularity-bound'
        result['scalar_bounds'] = {name: envelope(value, args.output_bits) for name, value in scalar.items()}
        for name, value in scalar.items():
            print(name, value.str(14), flush=True)
    except (ArithmeticError, ValueError, ZeroDivisionError) as error:
        result.update(verdict='inconclusive-exact-reference', diagnostic=str(error))
    result['elapsed_seconds'] = time.monotonic()-started
    write_packet(args.out, result)
    print(result['verdict'], 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(core='joint-gamma-6765-rho32-v2-schur.json',
                             cross='full-gamma-defect-cross.json',
                             rank='joint-gamma-6765-rho32-v2-rank-floor.json',
                             shared='full-gamma-correlated-sign.json').items():
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
