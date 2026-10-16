"""Pay the missing infinite residual tail using accepted cached inputs.

Reuse the complete finite dual defect, inverse profiles, same-direction
far residual energy and existing normalized symbol polynomial/remainder.
Only scalar norm bounds are evaluated; no scalar/moment supplier, matrix
action, solve, optimizer or larger Fourier cutoff is run.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from flint import arb, ctx, fmpq

from arb_arithmetic import Arithmetic
from arb_far_coupling import SEED_SHA256
from arb_inverse_residual_obstruction import HERE, M, N
from diagonal_floor_witness import interval_ball
from gamma_boundary_update import read_json
from joint_gamma_boundary import matrix_in
from structured_inverse_residual_interval import envelope


def consume(args):
    started = time.monotonic()
    finite, finite_hash = read_json(args.finite)
    defect, defect_hash = read_json(args.defect)
    direction, direction_hash = read_json(args.direction)
    pilot, pilot_hash = read_json(args.pilot)
    seed, seed_hash = read_json(args.seed)
    if (finite['schema'], defect['schema'], direction['schema'], seed['schema'],
            defect['metadata'], defect['finite_action_sha256'],
            finite['metadata']['direction_sha256'],
            finite['metadata']['input_sha256']['inverse-residual-pilot.json.zlib'],
            finite['dual_first'], finite['dual_last'], seed_hash,
            seed['c'], seed['M'], seed['first_mode'], seed['stirling_order']) != (
            'ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1', 'ACTUAL_C27_FINITE_DUAL_DEFECT_V1',
            'ACTUAL_C27_JOINT_GAMMA_DIRECTION_V1', 'ACTUAL_WEIL_NORMALIZED_TAIL_MOMENTS_V1',
            finite['metadata'], finite_hash, direction_hash, pilot_hash,
            2049, M, SEED_SHA256, 27, M, M+1, 32):
        raise ValueError('Require the same accepted finite profiles, current direction and tail seed.')
    p = matrix_in(defect['F_inverse_z'], defect['grid_bits'])
    q = matrix_in(defect['D_inverse_difference'], defect['grid_bits'])
    if (p.nrows(), p.ncols(), q.nrows(), q.ncols()) != (768, 1, 768, 1):
        raise ValueError('Incomplete accepted inverse profiles.')
    w = [arb(fmpq(int(x), 1 << finite['dual_bits'])) for x in finite['dual_integers']]
    if len(w) != M-2048:
        raise ValueError('Incomplete accepted exact dual vector.')
    theta = arb(fmpq(int(defect['theta_integer']), 1 << defect['theta_bits']))
    Dfar = arb(fmpq(int(pilot['far_weight_numerator']),
                    int(pilot['far_weight_denominator'])))
    arithmetic = Arithmetic(27, args.output_bits)
    d, pi, k = arithmetic.d, arithmetic.pi, 2*arithmetic.pi/arithmetic.d
    center = (arb(fmpq(1, 4))-1/arb(27).sqrt()/(1-arb(fmpq(1, 729)))
              +2*(arithmetic.H-1))
    coefficients = [interval_ball(x, seed['grid_bits']) for x in seed['normalized_polynomial']]
    if len(coefficients) != 64 or not coefficients[0].overlaps(-pi/4):
        raise ValueError('Wrong accepted smooth GP polynomial.')
    epsilon = interval_ball(seed['uniform_symbol_error'], seed['analytic_grid_bits']).upper()
    A, ratio = M+1, arb(fmpq(M, M+1))
    # C0=-pi/4 is exact analytically: its exported interval is NOT multiplied
    # by unbounded m. Higher terms and the source-proved decaying remainder
    # have maximal m-weighted absolute bounds at m=A.
    tail_rho = abs(-k*M*coefficients[1]-center)+k*A*(
        sum((abs(coefficients[j])*ratio**j for j in range(2, len(coefficients))), arb(0))
        +epsilon)
    gp = {row['n']: interval_ball(row['s'], finite['grid_bits'])
          for row in finite['boundary_scalars']}
    if set(gp) != set(range(N+1, M+1)) or len(gp) != len(finite['boundary_scalars']):
        raise ValueError('Incomplete accepted finite full-rate symbol.')
    rho = {n: abs(-k*n*gp[n]-pi*k*n/4-center) for n in gp}
    gamma = d/pi**2

    def gp_tail_norm(values, first):
        one = sum(map(abs, values), arb(0))
        remainder = arb(0)
        for i, x in enumerate(values):
            n = first+i
            a = M+1-n
            if a < 1:
                raise ValueError('Finite source must end before the far region.')
            square_sum = (arb(fmpq(1, a*a))+arb(fmpq(1, a)))/(M+1+n)**2
            remainder += abs(x)*(tail_rho+rho[n])*square_sum.sqrt()
        main = one/(2*arb(M).sqrt())
        return Dfar.sqrt()*(main+gamma*remainder), main, gamma*remainder

    wnorm, wmain, werror = gp_tail_norm(w, 2049)
    qnorm, qmain, qerror = gp_tail_norm(list(q.entries()), N+1)
    # Same old actual-symbol budget and separated divided-difference bound.
    budget = 2*arithmetic.H+sum((weight for weight, _ in arithmetic.primes), arb(0))
    difference_one = sum(map(abs, (q-p).entries()), arb(0))
    actual_norm = Dfar.sqrt()*2*budget/pi*difference_one/arb(M-2048).sqrt()
    residual_far = interval_ball(direction['actual_far_residual_energy_upper'], direction['grid_bits'])
    if not residual_far > 0 or not Dfar > 0:
        raise ValueError('Positive accepted far allowance and original weight required.')
    total_action_norm = wnorm+actual_norm+qnorm
    tail_norm_bound = residual_far.upper().sqrt()+abs(theta)*total_action_norm
    tail_debit = (tail_norm_bound**2).upper()
    finite_lower = interval_ball(defect['finite_optimum_lower'], defect['grid_bits'])
    lower = finite_lower-tail_debit
    fields = dict(smooth_tail_rho_bound=tail_rho, actual_symbol_budget=budget,
                  gp_w_tail_norm_upper=wnorm, gp_w_hankel_term=wmain,
                  gp_w_remainder_term=werror, gp_q_tail_norm_upper=qnorm,
                  gp_q_hankel_term=qmain, gp_q_remainder_term=qerror,
                  actual_mixed_tail_norm_upper=actual_norm,
                  effective_action_tail_norm_upper=total_action_norm,
                  residual_tail_norm_upper=tail_norm_bound,
                  additional_far_defect_upper=tail_debit,
                  finite_reference_lower=finite_lower,
                  comparison_trial_margin_lower=lower)
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_PRIMAL_TAIL_V1', metadata=finite['metadata'],
                  finite_action_sha256=finite_hash, finite_defect_sha256=defect_hash,
                  direction_sha256=direction_hash, pilot_sha256=pilot_hash,
                  tail_seed_sha256=seed_hash, grid_bits=args.output_bits,
                  working_bits=args.precision,
                  support_trial='theta*(D^-1(z-b)-F^-1z), with the same accepted exact theta and inverse profiles',
                  bounds={key: envelope(value, args.output_bits) for key, value in fields.items()},
                  verdict='same-direction-full-Gamma-comparison-positive' if lower > 0 else 'inconclusive',
                  conditions=defect['conditions']+[
                      'The accepted smooth GP polynomial has exact constant coefficient -pi/4 and the source-proved pointwise remainder decay.',
                      'The matching saved far residual energy applies to this exact retained v and selected inverse x0.',
                      'Finite support correction and inverse-residual square completion are valid on the inherited actual common form/domain.'],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started,
                  scope='Pay only the previously missing infinite residual tail for the explicit same-direction support trial. '
                    'A positive lower supplies this scalar sufficient comparison, not a complete retained-matrix sign, '
                    'a new generic theorem/algorithm, Lean verification, cofinal positivity or RH.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['verdict'], 'lower', lower.str(24), 'far debit', tail_debit.str(24),
          'seconds', time.monotonic()-started, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(finite='full-gamma-dual-32768-finite.json.zlib',
                             defect='full-gamma-dual-32768-defect.json.zlib',
                             direction='joint-gamma-6765-direction.json',
                             pilot='inverse-residual-pilot.json.zlib',
                             seed='tail-moments.json.zlib').items():
        parser.add_argument('--'+name, type=Path, default=HERE/default)
    parser.add_argument('--precision', type=int, default=1024)
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.output_bits != 640 or args.precision < 832:
        parser.error('Use a fresh output, the inherited640-bit grid and at least192 guard bits.')
    ctx.prec = args.precision
    consume(args)


if __name__ == '__main__':
    main()
