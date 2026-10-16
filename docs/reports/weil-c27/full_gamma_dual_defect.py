"""Apply existing finite inverse-defect bounds to the saved complete dual.

Only the missing residual and effective dual-action profiles are supplied.
Squared integer nodes turn the existing divided-difference action into
ordinary Toeplitz/Hankel polynomial products, evaluated by FLINT. This is
an application adapter, not a new convolution or inverse algorithm.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from flint import arb, arb_mat, arb_poly, ctx, fmpq

from arb_inverse_residual_obstruction import HERE, M, N
from arb_inverse_residual_pilot import trial_integer
from diagonal_floor_witness import interval_ball
from full_gamma_dual import boundary_entry, source_inputs
from full_gamma_dual_middle import finite_right_action
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope


class SquaredNodeAction:
    """Consume two existing polynomial products per squared-node K action."""

    def __init__(self, first, last):
        if not 0 < first <= last:
            raise ValueError('A consecutive positive finite node range is required.')
        self.first, self.count = first, last-first+1
        q = self.count
        self.toeplitz = arb_poly([arb(fmpq(1, k)) if k else arb(0)
                                 for k in range(1-q, q)])
        self.hankel = arb_poly([arb(fmpq(1, 2*first+k)) for k in range(2*q-1)])

    def kernel(self, values):
        q = self.count
        if len(values) != q:
            raise ValueError('Incomplete finite squared-node source.')
        toeplitz = self.toeplitz*arb_poly(values)
        hankel = self.hankel*arb_poly(list(reversed(values)))
        return [(toeplitz[q-1+i]+hankel[q-1+i]-values[i]/(2*(self.first+i)))
                /(2*(self.first+i)) for i in range(q)]

    def apply(self, values, phi, diagonal, gamma):
        q = self.count
        if len(phi) != q or len(diagonal) != q:
            raise ValueError('Incomplete finite symbol or independent diagonal.')
        u = [(-1)**(self.first+i)*values[i] for i in range(q)]
        ku = self.kernel(u)
        kphi = self.kernel([phi[i]*u[i] for i in range(q)])
        return [(-1)**(self.first+i)*(diagonal[i]*u[i]
                +gamma*(phi[i]*ku[i]-kphi[i])) for i in range(q)]


def produce(args):
    started = time.monotonic()
    (arithmetic, data, pilot, old, generators, selected, _, _, GP,
     F, vq, provenance) = source_inputs(args)
    saved, finite_hash = read_json(args.finite)
    bound, bound_hash = read_json(args.bound)
    if (saved['schema'], saved['metadata'], saved['dual_first'], saved['dual_last'],
            bound['schema'], bound['metadata'], bound['finite_action_sha256'],
            bound['D_floor']) != (
            'ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1', provenance, 2049, M,
            'ACTUAL_C27_FULL_GAMMA_DUAL_BOUND_V1', provenance, finite_hash, '1/632'):
        raise ValueError('Require the accepted same-direction complete middle action and bound.')
    if len(saved['dual_integers']) != M-2048:
        raise ValueError('Incomplete saved exact dual vector.')
    bits = saved['grid_bits']
    gp = {row['n']: (interval_ball(row['s'], bits), interval_ball(row['a'], bits))
          for row in saved['boundary_scalars']}
    periodic = {row['n']: interval_ball(row['periodic'], bits)
                for row in saved['boundary_scalars']}
    modes = list(range(N+1, M+1))
    if set(gp) != set(modes) or len(gp) != len(saved['boundary_scalars']):
        raise ValueError('Incomplete or repeated full-rate scalar nodes.')
    weights = [arb(fmpq(int(x), 1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    B0 = arb_mat([[boundary_entry(arithmetic, m, n, gp, periodic) for n in selected]
                  for m in selected])
    for i, m in enumerate(selected):
        B0[i, i] += 1/weights[m-N-1]
    D = F-B0
    z, b = matrix_in(saved['z'], bits), matrix_in(saved['b'], bits)
    if (z.nrows(), z.ncols(), b.nrows(), b.ncols()) != (768, 1, 768, 1):
        raise ValueError('Incomplete saved common mixed actions.')
    # These inverse vectors were not archived by the scalar bound consumer.
    # Its same certified solve is reused; no new inverse/floor theorem enters.
    fz = F.solve(z, algorithm='precond')
    dk = D.solve(z-b, algorithm='precond')
    w = [arb(fmpq(int(x), 1 << saved['dual_bits'])) for x in saved['dual_integers']]
    action = SquaredNodeAction(N+1, M)
    gamma = arithmetic.d/arithmetic.pi**2
    omega = [2*arithmetic.pi*n/arithmetic.d for n in modes]
    phi_gp = [-omega[i]*gp[n][0] for i, n in enumerate(modes)]
    diag_gp = [gp[n][1]-gp[n][0]/(arithmetic.pi*n)-periodic[n] for n in modes]
    phi_actual = [-omega[i]*data[n][0] for i, n in enumerate(modes)]
    diag_actual = [data[n][1]-data[n][0]/(arithmetic.pi*n) for n in modes]
    # M_J w = Lambda_J w + G_JJ w - H_JE F^-1 z
    #                       + (H_JE-G_JE) D^-1(z-b).
    # Both complete mixed terms are included by the following two actions.
    g_action = action.apply([-x for x in dk.entries()]+w, phi_gp, diag_gp, gamma)
    h_action = action.apply(list((dk-fz).entries())+[arb(0)]*len(w),
                            phi_actual, diag_actual, gamma)
    offset = len(selected)
    mw = [w[i]/weights[offset+i]+g_action[offset+i]+h_action[offset+i]
          for i in range(len(w))]
    GO = matrix_in(generators['GO'], generators['grid_bits'])
    kap = interval_ball(old['kappa'], old['grid_bits'])
    outside = list(range(2049, M+1))
    rmat = finite_right_action(outside, lambda m, j: GO[m-2049, j], GP,
                              arb_mat([[arb(x)] for x in vq]), list(range(N+1)),
                              N, 256, kap)
    r = list(rmat.entries())
    actual_a = sum((w[i]*r[i] for i in range(len(w))), arb(0))
    actual_h = sum((w[i]*mw[i] for i in range(len(w))), arb(0))
    if not actual_a.overlaps(interval_ball(saved['a'], bits)):
        raise ArithmeticError('Complete residual disagrees with accepted common scalar action.')
    if not actual_h.overlaps(interval_ball(bound['denominator'], bound['grid_bits'])):
        raise ArithmeticError('Effective dual action disagrees with accepted denominator.')
    if not actual_h > 0:
        raise ArithmeticError('Positive same-vector dual denominator is required.')
    theta_integer = trial_integer(actual_a/actual_h, 512)
    theta = arb(fmpq(theta_integer, 1 << 512))
    linear = 2*theta*actual_a-theta**2*actual_h
    error = sum((weights[offset+i]*(r[i]-theta*mw[i])**2
                 for i in range(len(w))), arb(0))
    Q = interval_ball(saved['Q'], bits)
    upper, lower = Q-linear, Q-linear-error
    verdict = ('complete-Gamma-comparison-excluded' if upper < 0 else
               'this-complete-finite-dual-space-insufficient' if lower > 0 else
               'finite-dual-optimization-defect-unresolved')
    # Direct common-object row checks use cached entries, not new suppliers.
    samples = [2049, 4096, M]
    for m in samples:
        i = m-2049
        direct_r = sum((kap*(GO[i, 0]*GP[n, 1]-GO[i, 1]*GP[n, 0])
                        /(m*m-n*n)*arb(vq[n]) for n in range(N+1)), arb(0))
        direct_m = w[i]/weights[offset+i]
        direct_m += sum((boundary_entry(arithmetic, m, n, gp, periodic)*w[j]
                         for j, n in enumerate(outside)), arb(0))
        direct_m += sum((arithmetic.even_entry(m, n, data)*(dk[j, 0]-fz[j, 0])
                          -boundary_entry(arithmetic, m, n, gp, periodic)*dk[j, 0]
                          for j, n in enumerate(selected)), arb(0))
        if not r[i].overlaps(direct_r) or not mw[i].overlaps(direct_m):
            raise ArithmeticError('Complete profile disagrees with an actual direct row.')
    result = dict(schema='ACTUAL_C27_FINITE_DUAL_DEFECT_V1', metadata=provenance,
                  finite_action_sha256=finite_hash, inherited_bound_sha256=bound_hash,
                  working_bits=args.precision, grid_bits=args.output_bits,
                  dual_first=2049, dual_last=M,
                  theta_bits=512, theta_integer=str(theta_integer),
                  Q=saved['Q'], a=envelope(actual_a, args.output_bits),
                  denominator=envelope(actual_h, args.output_bits),
                  dual_trial_gain=envelope(linear, args.output_bits),
                  weighted_residual_defect=envelope(error, args.output_bits),
                  finite_optimum_lower=envelope(lower, args.output_bits),
                  finite_optimum_upper=envelope(upper, args.output_bits),
                  residual_profile=matrix_out(rmat, args.output_bits),
                  effective_dual_action=matrix_out(arb_mat([[x] for x in mw]), args.output_bits),
                  F_inverse_z=matrix_out(fz, args.output_bits),
                  D_inverse_difference=matrix_out(dk, args.output_bits),
                  direct_sample_modes=samples, verdict=verdict,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started,
                  conditions=bound['conditions']+[
                      'D>0 and T>=Lambda imply M_J>=Lambda_J by the finite block factorization.',
                      'The inherited floor and exact dyadic w, Q, z, b, scalar array and original W are sound.'],
                  scope='The inverse-defect enclosure concerns the entire current finite dual space '
                    '2049..32768 for the same full-Gamma comparison. A positive lower only rules '
                    'out this finite restriction; no actual Schur, larger-space, cofinal or RH sign. '
                    'No new general algorithm, numerical supplier or Lean claim.')
    write_packet(args.out, result)
    print(verdict, 'lower', lower.str(24), 'upper', upper.str(24),
          'defect', error.str(24), 'seconds', time.monotonic()-started, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(target='structured-joint-2048-target.json.zlib',
                             pivots='structured-joint-2048-pivots.json',
                             support='support-2048-finite-action.json.zlib',
                             generators='joint-residual-generators-2048.json.zlib',
                             direction='joint-gamma-6765-direction.json',
                             finite='full-gamma-dual-32768-finite.json.zlib',
                             bound='full-gamma-dual-32768-bound.json').items():
        p.add_argument('--'+name, type=Path, default=HERE/default)
    p.add_argument('--precision', type=int, default=1024)
    p.add_argument('--output-bits', type=int, default=640)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if args.out.exists() or args.output_bits != 640 or args.precision < 832:
        p.error('Use a fresh output, the inherited640-bit grid and at least192 guard bits.')
    ctx.prec = args.precision
    produce(args)


if __name__ == '__main__':
    main()
