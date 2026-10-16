"""Extend the accepted finite dual vector to every archived middle row.

Reuse its exact prefix, Q and common actions. Existing finite Cauchy
products supply only the new tail actions. A finite-node triangle
envelope around a Hankel kernel makes the complete boundary quadratic
available through FLINT polynomial multiplication. No asymptotic or
infinite-grid operator bound is assumed.
"""
import argparse
import hashlib
from pathlib import Path
import time

from flint import arb, arb_mat, arb_poly, ctx, fmpq

from arb_inverse_residual_obstruction import HERE, M, N
from arb_inverse_residual_pilot import trial_integer
from diagonal_floor_witness import interval_ball
from full_gamma_dual import full_boundary_scalars, proposal_channels, source_inputs
from gamma_boundary_update import read_json
from gamma_channel_resolution import cauchy_projection
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope


def finite_right_action(modes, row_generator, column_generator, coefficients,
                        columns, scale, order, kappa):
    """Apply the existing finite disjoint geometric Cauchy expansion."""
    a = [column_generator[i, 1]*coefficients[i, 0] for i in range(len(columns))]
    b = [column_generator[i, 0]*coefficients[i, 0] for i in range(len(columns))]
    absolute_a, absolute_b = sum(map(abs, a), arb(0)), sum(map(abs, b), arb(0))
    power = [arb(1)]*len(columns)
    ratios = [arb(fmpq(n*n, scale*scale)) for n in columns]
    ma, mb = [], []
    for _ in range(order):
        ma.append(sum((a[i]*power[i] for i in range(len(columns))), arb(0)))
        mb.append(sum((b[i]*power[i] for i in range(len(columns))), arb(0)))
        power = [power[i]*ratios[i] for i in range(len(columns))]
    pa, pb = arb_poly(ma), arb_poly(mb)
    result = arb_mat(len(modes), 1)
    for i, m in enumerate(modes):
        t = arb(fmpq(scale*scale, m*m))
        gm0, gm1 = row_generator(m, 0), row_generator(m, 1)
        value = kappa/(m*m)*(gm0*pa(t)-gm1*pb(t))
        error = kappa/(m*m)*t**order/(1-t)*(abs(gm0)*absolute_a+abs(gm1)*absolute_b)
        result[i, 0] = value+arb(0, error.upper())
    return result


def hankel_quadratic(w, first):
    # Existing FLINT polynomial product; all coefficients are real here.
    p = arb_poly([arb(0)]*first+[(-1)**(first+i)*w[i, 0] for i in range(w.nrows())])
    squared = p*p
    return sum((squared[k]/(2*k) for k in range(2*first, len(squared))), arb(0))


def produce(args):
    started = time.monotonic()
    (arithmetic, data, pilot, old, generators, selected, GE, _, GP, _, vq,
     provenance) = source_inputs(args)
    prefix, prefix_hash = read_json(args.prefix)
    trial, trial_hash = read_json(args.proposal)
    boundary, boundary_hash = read_json(args.boundary)
    refinement, refinement_hash = read_json(args.refinement)
    if (prefix['schema'], prefix['metadata'], prefix['dual_first'], prefix['dual_last'],
            prefix['proposal_sha256'], prefix['boundary_sha256'], prefix['refinement_sha256']) != (
            'ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1', provenance, 2049, 4096,
            trial_hash, boundary_hash, refinement_hash):
        raise ValueError('Require the accepted same-direction2049..4096 dual prefix.')
    bits = prefix['grid_bits']
    tail = list(range(prefix['dual_last']+1, M+1))
    v = arb_mat([[arb(x)] for x in vq])
    t = arb_mat([[arb(fmpq(int(x), 1 << trial['trial_correction_bits']))]
                 for x in trial['trial_correction_integers']])
    zc = arb_mat([[arb(fmpq(int(x), 1 << trial['channel_trial_bits']))]
                  for x in trial['channel_trial_integers']])
    GO = matrix_in(generators['GO'], generators['grid_bits'])
    kap = interval_ball(old['kappa'], old['grid_bits'])
    residual = finite_right_action(tail, lambda m, j: GO[m-2049, j], GP, v,
                                   list(range(N+1)), N, 256, kap)
    original = lambda m, j: arb((-1)**m)*(m*data[m][0] if j else 1)
    zt = finite_right_action(tail, original, GE, t, selected, 2048, 384, kap)
    channels = proposal_channels(arithmetic, boundary['bins'], tail)
    proposal = residual-zt-channels*zc
    weights = [arb(fmpq(int(x), 1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    shift, dual_bits = prefix['proposal_power_two_shift'], prefix['dual_bits']
    tail_integers = [trial_integer(weights[m-N-1]*proposal[i, 0]*arb(2)**shift, dual_bits)
                     for i, m in enumerate(tail)]
    wt = arb_mat([[arb(fmpq(x, 1 << dual_bits))] for x in tail_integers])
    wp = arb_mat([[arb(fmpq(int(x), 1 << dual_bits))] for x in prefix['dual_integers']])
    w = arb_mat([[x] for x in wp.entries()+wt.entries()])
    a = interval_ball(prefix['a'], bits)+(wt.transpose()*residual)[0, 0]
    weighted = arb_mat(M-N, 1)
    for i, m in enumerate(tail):
        weighted[m-N-1, 0] = wt[i, 0]
    ztail, z_geometry = cauchy_projection(weighted, original, GE, tail, selected, kap, 2048, 384)
    z = matrix_in(prefix['z'], bits)+ztail.transpose()
    gp, periodic = full_boundary_scalars(arithmetic, tail, data)
    gp.update({row['n']: (interval_ball(row['s'], bits), interval_ball(row['a'], bits))
               for row in prefix['boundary_scalars']})
    periodic.update({row['n']: interval_ball(row['periodic'], bits) for row in prefix['boundary_scalars']})
    prime_generator = lambda m, j: arb((-1)**m)*(m*(data[m][0]-gp[m][0]) if j else 1)
    Gprime = arb_mat([[prime_generator(n, 0), prime_generator(n, 1)] for n in selected])
    utail, u_geometry = cauchy_projection(weighted, prime_generator, Gprime,
                                         tail, selected, kap, 2048, 384)
    u = matrix_in(prefix['z'], bits)-matrix_in(prefix['b'], bits)+utail.transpose()
    b = z-u
    print('new full-middle common thin actions completed; seconds', time.monotonic()-started, flush=True)
    # The constant only centers a finite array; its value is not a lemma.
    constant = (arb(fmpq(1, 4))-1/arb(27).sqrt()/(1-arb(fmpq(1, 729)))
                +2*(arithmetic.H-1))
    emax, diagonal_max = arb(0), arb(0)
    first = prefix['dual_first']
    for n in range(first, M+1):
        omega = 2*arithmetic.pi*n/arithmetic.d
        error = -omega*gp[n][0]-arithmetic.pi*omega/4-constant
        emax = max(emax, abs(error).upper())
        if n > prefix['dual_last']:
            diagonal_error = gp[n][1]-gp[n][0]/(arithmetic.pi*n)-periodic[n]-arb(fmpq(1, 4*n))
            diagonal_max = max(diagonal_max, abs(diagonal_error).upper())
    norm1 = sum(map(abs, w.entries()), arb(0))
    norm2 = sum((x*x for x in w.entries()), arb(0))
    pnorm1 = sum(map(abs, wp.entries()), arb(0))
    pnorm2 = sum((x*x for x in wp.entries()), arb(0))
    entry_error = 8/arithmetic.d*emax/((2*arithmetic.pi/arithmetic.d)**2*(2*first+1))
    offdiagonal_mass = norm1**2-norm2-pnorm1**2+pnorm2
    error = entry_error*offdiagonal_mass+diagonal_max*(norm2-pnorm2)
    if not error >= 0:
        raise ArithmeticError('Nonnegative finite Hankel error required.')
    diagonal_tail = sum((wt[i, 0]**2/weights[m-N-1] for i, m in enumerate(tail)), arb(0))
    hk_full, hk_prefix = hankel_quadratic(w, first), hankel_quadratic(wp, first)
    center = interval_ball(prefix['beta'], bits)+diagonal_tail+hk_full-hk_prefix
    beta = center+arb(0, error.upper())
    # Existing exact forms also check the finite geometric actions at endpoints.
    sample_modes = [tail[0], 16384, M]
    for m in sample_modes:
        i = m-tail[0]
        direct_r = sum((kap*(GO[m-2049, 0]*GP[n, 1]-GO[m-2049, 1]*GP[n, 0])
                        /(m*m-n*n)*v[n, 0] for n in range(N+1)), arb(0))
        direct_zt = sum((arithmetic.even_entry(m, n, data)*t[j, 0]
                         for j, n in enumerate(selected)), arb(0))
        if not residual[i, 0].overlaps(direct_r) or not zt[i, 0].overlaps(direct_zt):
            raise ArithmeticError('New finite tail expansion disagrees with actual direct row.')
    packet = dict(schema=prefix['schema'], metadata=provenance,
                  parent_dual_action_sha256=prefix_hash, proposal_sha256=trial_hash,
                  boundary_sha256=boundary_hash, refinement_sha256=refinement_hash,
                  grid_bits=args.output_bits, working_bits=args.precision,
                  dual_first=first, dual_last=M, dual_bits=dual_bits,
                  dual_integers=prefix['dual_integers']+list(map(str, tail_integers)),
                  proposal_power_two_shift=shift, Q=prefix['Q'],
                  a=envelope(a, args.output_bits), z=matrix_out(z, args.output_bits),
                  b=matrix_out(b, args.output_bits), beta=envelope(beta, args.output_bits),
                  boundary_scalars=prefix['boundary_scalars']+[
                      dict(n=n, s=envelope(gp[n][0], args.output_bits),
                           a=envelope(gp[n][1], args.output_bits),
                           periodic=envelope(periodic[n], args.output_bits)) for n in tail],
                  hankel=dict(center=envelope(center, args.output_bits),
                              error_bound=envelope(error, args.output_bits),
                              scalar_deviation=envelope(emax, args.output_bits),
                              diagonal_deviation=envelope(diagonal_max, args.output_bits),
                              offdiagonal_mass=envelope(offdiagonal_mass, args.output_bits),
                              new_diagonal_mass=envelope(norm2-pnorm2, args.output_bits),
                              finite_node_first=first, finite_node_last=M,
                              method='Finite-node triangle envelope; no asymptotic or infinite-grid bound.'),
                  new_residual_action_order=256, new_support_action_order=384,
                  direct_sample_modes=sample_modes, z_geometry=z_geometry, u_geometry=u_geometry,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started,
                  scope='Extend the same exact finite dual prefix through every archived middle mode. '
                    'Reuse current Q and prefix actions; existing Cauchy products and FLINT convolution only. '
                    'Full-rate boundary quadratic has a paid finite-node error; no discarded Fourier residual. '
                    'Requires inherited actual form/domain and sound inputs; no actual Schur sign, Lean, cofinal positivity or RH result.')
    write_packet(args.out, packet)
    print('full-middle dual archived; beta error', error.str(20),
          'seconds', time.monotonic()-started, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name, default in dict(target='structured-joint-2048-target.json.zlib',
                             pivots='structured-joint-2048-pivots.json',
                             support='support-2048-finite-action.json.zlib',
                             generators='joint-residual-generators-2048.json.zlib',
                             direction='joint-gamma-6765-direction.json',
                             proposal='gamma-support-pole-contrasts-depth1-optimum.json',
                             boundary='joint-gamma-6765-finite.json.zlib',
                             refinement='gamma-rate-contrasts-depth1.json.zlib',
                             prefix='full-gamma-dual-4096-finite.json.zlib').items():
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
