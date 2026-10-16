"""Restrict the full Gamma/pole inverse variational formula to [E,w].

Reuse actual Fourier scalars and the saved selected inverse/residual data.
Only the new complete-direction quadratic and finite common dual actions
are produced. This is an application of finite block elimination, not a
new inverse algorithm or an actual Weil-sign certificate.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

from flint import acb, arb, arb_mat, ctx, fmpq

from arb_arithmetic import Arithmetic
from arb_inverse_residual_obstruction import HERE, M, N, inputs
from arb_inverse_residual_pilot import trial_integer
from diagonal_floor_witness import interval_ball
from gamma_boundary_update import read_json
from gamma_channel_resolution import refinement_tree
from joint_gamma_boundary import matrix_in, matrix_out, original_target, write_packet
from structured_inverse_residual_interval import eliminate, envelope


def source_inputs(args):
    old, target_hash, _ = original_target(args)
    support, support_hash = read_json(args.support)
    generators, generator_hash = read_json(args.generators)
    candidate, candidate_hash = read_json(args.direction)
    if (support['schema'], generators['schema'], candidate['schema'],
            support['metadata']['target_packet_sha256'],
            generators['metadata']['target_packet_sha256'],
            generators['support_action_sha256'], candidate['original_target_sha256']) != (
            'ACTUAL_C27_SUPPORT_ACTION_V1', 'ACTUAL_C27_RESIDUAL_GENERATOR_V1',
            'ACTUAL_C27_JOINT_GAMMA_DIRECTION_V1', target_hash, target_hash,
            support_hash, target_hash):
        raise ValueError('Require the same complete target, selected block and direction.')
    scalar, rows, pilot, _ = inputs()
    arithmetic = Arithmetic(27, scalar['grid_bits'])
    selected = list(range(N+1, 2049))
    if old['metadata']['last_mode'] != selected[-1]:
        raise ValueError('Require the existing selected support through2048.')
    vq = [fmpq(int(x), 1 << candidate['direction_dyadic_bits'])
          for x in candidate['direction_dyadic_integers']]
    if len(vq) != N+1:
        raise ValueError('Incomplete retained direction.')
    data = {n: (interval_ball(row['s'], scalar['grid_bits']),
                interval_ball(row['a'], scalar['grid_bits'])) for n, row in rows.items()}
    GE = matrix_in(support['GE'], support['grid_bits'])
    inverse_GE = matrix_in(generators['selected_inverse_generator'], generators['grid_bits'])
    GP = arb_mat([[interval_ball(x, old['grid_bits']) for x in row[:2]]
                  for row in old['generators']])
    if (GE.nrows(), GE.ncols(), inverse_GE.nrows(), inverse_GE.ncols(),
            GP.nrows(), GP.ncols()) != (768, 2, 768, 2, N+1, 2):
        raise ValueError('Incomplete exact selected inverse interface.')
    F = arb_mat([[arithmetic.even_entry(m, n, data) for n in selected] for m in selected])
    for i in range(len(selected)):
        if not F[i, i].overlaps(interval_ball(support['F_diagonal'][i], support['grid_bits'])):
            raise ValueError('Selected diagonal changed.')
    provenance = dict(target_sha256=target_hash, support_sha256=support_hash,
                      generators_sha256=generator_hash, direction_sha256=candidate_hash,
                      input_sha256=old['metadata']['input_sha256'], N=N, M=M,
                      last_mode=selected[-1], basis=old['metadata']['basis'])
    return arithmetic, data, pilot, old, generators, selected, GE, inverse_GE, GP, F, vq, provenance


def full_boundary_scalars(arithmetic, modes, data):
    """Remove prime terms and the periodic Gamma diagonal from cached scalars.

    The independent even diagonal includes -sGP/(pi*n). No discrete
    differentiation or finite-channel trace bound enters this interface.
    All modes here are positive, in the same unphased orthonormal basis.
    """
    gp, periodic = {}, {}
    for n in modes:
        omega = 2*arithmetic.pi*n/arithmetic.d
        sine, diagonal = arb(0), arb(0)
        for weight, shift in arithmetic.primes:
            sine += weight*(omega*shift).sin()
            diagonal += 2*weight*(1-shift/arithmetic.d)*(omega*shift).cos()
        gp[n] = (data[n][0]+sine, data[n][1]+diagonal)
        periodic[n] = acb(arb(fmpq(1, 4)), omega/2).digamma().real-arithmetic.pi.log()
    return gp, periodic


def boundary_entry(arithmetic, m, n, gp, periodic):
    value = arithmetic.even_entry(m, n, gp)
    return value-periodic[n] if m == n else value


def proposal_channels(arithmetic, bins, modes):
    """The saved35 channel shapes, evaluated only on the new dual rows.

    These are the same mean/contrast formulas; their sole purpose here is
    to propose w. The certificate uses the archived exact dyadic w only.
    """
    _, details = refinement_tree(bins, 1)
    endpoints = {x for a, b in bins for x in (a, b)} | {d['split'] for d in details}
    coefficient = lambda a: 1-(-arb(fmpq(4*a+1, 2))*arithmetic.d).exp()
    U = arb_mat(len(modes), len(bins)+len(details))
    for i, n in enumerate(modes):
        im = arithmetic.pi*n/arithmetic.d
        values = {e: acb(arb(e)+arb(fmpq(1, 4)), im).digamma().real for e in endpoints}
        phase = (-1)**n
        for k, (a, b) in enumerate(bins):
            U[i, k] = phase*(coefficient(a)/((b-a)*arithmetic.d)).sqrt()*(values[b]-values[a])
        for k, detail in enumerate(details):
            a, mid, b = detail['left'], detail['split'], detail['right']
            p, q, length = mid-a, b-mid, b-a
            root = bins[detail['root']][0]
            U[i, len(bins)+k] = phase*(coefficient(root)/arithmetic.d).sqrt()*(
                arb(fmpq(q, p*length)).sqrt()*(values[mid]-values[a])
                -arb(fmpq(p, q*length)).sqrt()*(values[b]-values[mid]))
        U[i, 0] *= (arithmetic.d/4).exp()
    return U


def produce(args):
    started = time.monotonic()
    (arithmetic, data, pilot, old, generators, selected, GE, inverse_GE,
     GP, F, vq, provenance) = source_inputs(args)
    if not selected[-1] < args.dual_last <= M:
        raise ValueError('New dual rows must follow E and lie in the archived middle.')
    outside = list(range(selected[-1]+1, args.dual_last+1))
    v = arb_mat([[arb(x)] for x in vq])
    kap = interval_ball(old['kappa'], old['grid_bits'])
    # This is the cached two-generator inverse identity, not a new F solve.
    y = arb_mat([[sum((kap*(inverse_GE[i, 0]*GP[n, 1]-inverse_GE[i, 1]*GP[n, 0])
                       /(m*m-n*n)*v[n, 0] for n in range(N+1)), arb(0))]
                 for i, m in enumerate(selected)])
    # New quadratic for this v; the old support packet's qv belongs to another v.
    qa = arb(0)
    for m in range(N+1):
        qa += v[m, 0]**2*arithmetic.even_entry(m, m, data)
        qa += 2*v[m, 0]*sum((arithmetic.even_entry(m, n, data)*v[n, 0]
                             for n in range(m)), arb(0))
    Q = qa-(y.transpose()*F*y)[0, 0]
    print('new complete-direction Q', Q.str(20), 'seconds', time.monotonic()-started, flush=True)
    GO = matrix_in(generators['GO'], generators['grid_bits'])
    if (GO.nrows(), GO.ncols()) != (M-selected[-1], 2):
        raise ValueError('Incomplete saved residual generator.')
    r = arb_mat([[sum((kap*(GO[m-selected[-1]-1, 0]*GP[n, 1]
                            -GO[m-selected[-1]-1, 1]*GP[n, 0])/(m*m-n*n)*v[n, 0]
                       for n in range(N+1)), arb(0))] for m in outside])
    trial, trial_hash = read_json(args.proposal)
    boundary, boundary_hash = read_json(args.boundary)
    refinement, refinement_hash = read_json(args.refinement)
    if (trial['schema'], trial['direction_sha256'], trial['metadata']['target_sha256'],
            trial['metadata']['generators_sha256'], trial['retained_even_pole'],
            trial['channel_count'], trial['metadata']['boundary_sha256'],
            trial['metadata']['channel_refinement_sha256']) != (
            'ACTUAL_C27_GAMMA_SUPPORT_OPTIMUM_V1', provenance['direction_sha256'],
            provenance['target_sha256'], provenance['generators_sha256'], True, 35,
            boundary_hash, refinement_hash):
        raise ValueError('Proposal must be the same archived35-channel/pole trial.')
    t = arb_mat([[arb(fmpq(int(x), 1 << trial['trial_correction_bits']))]
                 for x in trial['trial_correction_integers']])
    zc = arb_mat([[arb(fmpq(int(x), 1 << trial['channel_trial_bits']))]
                  for x in trial['channel_trial_integers']])
    U = proposal_channels(arithmetic, boundary['bins'], outside)
    Z = arb_mat([[arithmetic.even_entry(m, n, data) for n in selected] for m in outside])
    proposal = r-Z*t-U*zc
    weights = [arb(fmpq(int(x), 1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    proposal = arb_mat([[weights[m-N-1]*proposal[i, 0]] for i, m in enumerate(outside)])
    amplitude = max(abs(float(x.mid())) for x in proposal.entries())
    if not math.isfinite(amplitude) or amplitude == 0:
        raise ArithmeticError('Finite nonzero proposal required.')
    shift = -math.floor(math.log2(amplitude))
    wi = [trial_integer(x*arb(2)**shift, args.dual_bits) for x in proposal.entries()]
    if not any(wi):
        raise ArithmeticError('Dual rounding erased the direction.')
    w = arb_mat([[arb(fmpq(x, 1 << args.dual_bits))] for x in wi])
    a, z = (w.transpose()*r)[0, 0], Z.transpose()*w
    modes = selected+outside
    gp, periodic = full_boundary_scalars(arithmetic, modes, data)
    B0w = arb_mat([[sum((boundary_entry(arithmetic, n, m, gp, periodic)*w[i, 0]
                         for i, m in enumerate(outside)), arb(0))] for n in selected])
    beta = arb(0)
    for i, m in enumerate(outside):
        beta += w[i, 0]**2*(boundary_entry(arithmetic, m, m, gp, periodic)+1/weights[m-N-1])
        beta += 2*w[i, 0]*sum((boundary_entry(arithmetic, m, n, gp, periodic)*w[j, 0]
                               for j, n in enumerate(outside[:i])), arb(0))
    packet = dict(schema='ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1', metadata=provenance,
                  proposal_sha256=trial_hash, boundary_sha256=boundary_hash,
                  refinement_sha256=refinement_hash, grid_bits=args.output_bits,
                  working_bits=args.precision, dual_first=outside[0], dual_last=outside[-1],
                  dual_bits=args.dual_bits, dual_integers=list(map(str, wi)),
                  proposal_power_two_shift=shift,
                  Q=envelope(Q, args.output_bits), a=envelope(a, args.output_bits),
                  z=matrix_out(z, args.output_bits), b=matrix_out(B0w, args.output_bits),
                  beta=envelope(beta, args.output_bits),
                  boundary_scalars=[dict(n=n, s=envelope(gp[n][0], args.output_bits),
                       a=envelope(gp[n][1], args.output_bits),
                       periodic=envelope(periodic[n], args.output_bits)) for n in modes],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started,
                  scope='New common finite dual actions for all Gamma rates and the even pole. '
                  'Cached actual scalars, selected inverse and residual generators reused; no Fourier tail debit '
                  'because only dual trials are restricted. No actual Schur sign, cofinal positivity, Lean or RH result.')
    digest = write_packet(args.out, packet)
    print('finite full-Gamma dual actions archived', digest,
          'seconds', time.monotonic()-started, flush=True)


def consume(args):
    started = time.monotonic()
    (arithmetic, data, pilot, old, _, selected, GE, _, _, F, _, provenance) = source_inputs(args)
    saved, packet_hash = read_json(args.finite)
    if saved['schema'] != 'ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1' or saved['metadata'] != provenance:
        raise ValueError('Wrong full-Gamma dual action inputs.')
    bits = saved['grid_bits']
    if len(saved['dual_integers']) != saved['dual_last']-saved['dual_first']+1 or not any(
            int(x) for x in saved['dual_integers']):
        raise ValueError('Incomplete or zero exact dual vector.')
    if saved['dual_first'] != selected[-1]+1 or saved['dual_last'] > M:
        raise ValueError('Dual vector must be disjoint from E.')
    gp = {row['n']: (interval_ball(row['s'], bits), interval_ball(row['a'], bits))
          for row in saved['boundary_scalars']}
    periodic = {row['n']: interval_ball(row['periodic'], bits) for row in saved['boundary_scalars']}
    if set(gp) != set(range(N+1, saved['dual_last']+1)):
        raise ValueError('Incomplete finite full-boundary scalar supply.')
    weights = [arb(fmpq(int(x), 1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    B0 = arb_mat([[boundary_entry(arithmetic, m, n, gp, periodic) for n in selected]
                  for m in selected])
    for i, m in enumerate(selected):
        B0[i, i] += 1/weights[m-N-1]
    D = F-B0
    # Separate diagonals, same rank-two displacement LDL, new matrix D.
    gd = arb_mat([[GE[i, 0], GE[i, 1]-(-1)**m*m*gp[m][0], arb(0), arb(0)]
                  for i, m in enumerate(selected)])
    slack = fmpq(1, 632)
    verdict, failure, pivots, order = eliminate(
        [D[i, i]-arb(slack) for i in range(len(selected))], gd, arb_mat(0, 0),
        interval_ball(old['kappa'], old['grid_bits']), [m*m for m in selected],
        pivoting='max-diagonal')
    if verdict != 'positive-pivots':
        raise ArithmeticError(f'New D floor unresolved: {verdict} at {failure}.')
    z, b = matrix_in(saved['z'], bits), matrix_in(saved['b'], bits)
    if (z.nrows(), z.ncols(), b.nrows(), b.ncols()) != (768, 1, 768, 1):
        raise ValueError('Incomplete common dual cross actions.')
    Q, a, beta = [interval_ball(saved[key], bits) for key in ('Q', 'a', 'beta')]
    fz = F.solve(z, algorithm='precond')
    dzb = D.solve(z-b, algorithm='precond')
    denominator = beta-(z.transpose()*fz)[0, 0]+((z-b).transpose()*dzb)[0, 0]
    if not denominator > 0:
        raise ArithmeticError('Dual optimum denominator must be strictly positive.')
    upper = Q-a*a/denominator
    result = dict(schema='ACTUAL_C27_FULL_GAMMA_DUAL_BOUND_V1', metadata=provenance,
                  finite_action_sha256=packet_hash, working_bits=args.precision,
                  grid_bits=args.output_bits, dual_first=saved['dual_first'],
                  dual_last=saved['dual_last'], D_floor=str(slack),
                  D_floor_pivot_modes=[selected[i] for i in order],
                  D_floor_pivots=[envelope(p, args.output_bits) for p in pivots],
                  Q=envelope(Q, args.output_bits), a=envelope(a, args.output_bits),
                  denominator=envelope(denominator, args.output_bits),
                  dual_subtraction=envelope(a*a/denominator, args.output_bits),
                  all_channel_all_trial_upper=envelope(upper, args.output_bits),
                  verdict='all-Gamma-pole-refinements-on-this-support-excluded' if upper < 0 else 'inconclusive',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started,
                  conditions=['Same actual form/domain, exterior, original orthonormal mass and saved complete direction.',
                    'Actual H>=T=Lambda+full Gamma boundary+even pole>0 with the same archived Lambda=W^-1.',
                    'Sound saved exact-inverse/residual and common finite directed enclosures.'],
                  scope='Finite common dual restriction gives an upper bound for every t on E and every '
                    'channel scheme BB*<=Gamma_boundary+even_pole with fixed Lambda and retained v. '
                    'Negative upper excludes this entire comparison family, not the actual Schur form. '
                    'Nonnegative upper is inconclusive. No Lean, cofinal positivity or RH result.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(result['verdict'], 'upper', upper.str(24), 'Q', Q.str(24),
          'denominator', denominator.str(24), 'seconds', time.monotonic()-started, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--produce', action='store_true')
    p.add_argument('--target', type=Path, default=HERE/'structured-joint-2048-target.json.zlib')
    p.add_argument('--pivots', type=Path, default=HERE/'structured-joint-2048-pivots.json')
    p.add_argument('--support', type=Path, default=HERE/'support-2048-finite-action.json.zlib')
    p.add_argument('--generators', type=Path, default=HERE/'joint-residual-generators-2048.json.zlib')
    p.add_argument('--direction', type=Path, default=HERE/'joint-gamma-6765-direction.json')
    p.add_argument('--proposal', type=Path, default=HERE/'gamma-support-pole-contrasts-depth1-optimum.json')
    p.add_argument('--boundary', type=Path, default=HERE/'joint-gamma-6765-finite.json.zlib')
    p.add_argument('--refinement', type=Path, default=HERE/'gamma-rate-contrasts-depth1.json.zlib')
    p.add_argument('--finite', type=Path)
    p.add_argument('--dual-last', type=int, default=4096)
    p.add_argument('--dual-bits', type=int, default=384)
    p.add_argument('--precision', type=int, default=1024)
    p.add_argument('--output-bits', type=int, default=640)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if args.out.exists() or args.precision < max(640, args.output_bits, args.dual_bits)+192:
        p.error('Use a fresh output path and at least192 working guard bits.')
    if not args.produce and args.finite is None:
        p.error('A cached finite dual action is required for the consumer.')
    ctx.prec = args.precision
    produce(args) if args.produce else consume(args)


if __name__ == '__main__':
    main()
