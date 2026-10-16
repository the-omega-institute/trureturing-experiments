"""Optimize the changed Gamma inverse comparison on a saved negative direction.

Reuse the old support Gram and mixed displacement identity. Produce only
missing Gamma/support thin projections and two mixed inverse generators.
Cached consumers evaluate a changed trial; they repeat no old producer.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

from flint import arb, arb_mat, arb_poly, ctx, fmpq

from arb_arithmetic import Arithmetic
from arb_far_coupling import reconstruct_factors
from arb_inverse_residual_obstruction import HERE, M, N, inputs
from arb_inverse_residual_pilot import trial_integer
from diagonal_floor_witness import interval_ball
from gamma_boundary_update import boundary_columns, read_json
from joint_gamma_boundary import matrix_in, matrix_out, original_target, write_packet
from structured_inverse_residual_interval import eliminate, envelope
from structured_support_test import gram_from_generators


def setup(args):
    old, th, metadata_bits = original_target(args)
    support, sh = read_json(args.support)
    generators, gh = read_json(args.generators)
    boundary, bh = read_json(args.boundary)
    last = old['metadata']['last_mode']
    if (support['schema'], generators['schema'], boundary['schema'], last,
            support['metadata']['target_packet_sha256'],
            generators['metadata']['target_packet_sha256'], boundary['target_packet_sha256'],
            generators['support_action_sha256'], boundary['generator_packet_sha256']) != (
            'ACTUAL_C27_SUPPORT_ACTION_V1', 'ACTUAL_C27_RESIDUAL_GENERATOR_V1',
            'ACTUAL_C27_JOINT_GAMMA_ACTION_V1', 2048, th, th, th, sh, gh):
        raise ValueError('Require the same complete target, support and Gamma action.')
    if boundary['input_sha256'] != old['metadata']['input_sha256']:
        raise ValueError('Gamma and original actual inputs differ.')
    selected = list(range(N+1, last+1))
    bits = support['grid_bits']
    GE = matrix_in(support['GE'], bits)
    TZ = matrix_in(support['T'], bits)
    diagonal = [interval_ball(x, bits) for x in support['gram_diagonal']]
    kap = interval_ball(old['kappa'], old['grid_bits'])
    F = arb_mat(len(selected), len(selected))
    for i, m in enumerate(selected):
        F[i, i] = interval_ball(support['F_diagonal'][i], bits)
        for j, n in enumerate(selected[:i]):
            F[i, j] = kap*(GE[i, 0]*GE[j, 1]-GE[i, 1]*GE[j, 0])/(m*m-n*n)
            F[j, i] = F[i, j]
    gram = gram_from_generators(GE, TZ, diagonal, kap, selected)
    inverse_GE = matrix_in(generators['selected_inverse_generator'], generators['grid_bits'])
    if (GE.nrows(), GE.ncols(), inverse_GE.nrows(), inverse_GE.ncols()) != (768, 2, 768, 2):
        raise ValueError('Incomplete selected generator action.')
    provenance = dict(target_sha256=th, support_sha256=sh, generators_sha256=gh,
                      boundary_sha256=bh, input_sha256=old['metadata']['input_sha256'],
                      N=N, M=M, last_mode=last, basis=old['metadata']['basis'])
    return old, metadata_bits, support, generators, boundary, selected, GE, TZ, F, gram, inverse_GE, provenance


def produce(args):
    started = time.monotonic()
    old, _, _, _, boundary, selected, GE, TZ, F, gram, inverse_GE, provenance = setup(args)
    # This is a new two-RHS action; the old F^-1 GE is consumed, not recomputed.
    inverse_mixed = F.solve(TZ-gram*inverse_GE, algorithm='precond')
    scalar, rows, pilot, _ = inputs()
    arithmetic = Arithmetic(27, scalar['grid_bits'])
    weights = [arb(fmpq(int(x), 1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    debit = arb(fmpq(int(pilot['far_weight_numerator']), int(pilot['far_weight_denominator'])))
    U, eta = boundary_columns(arithmetic, boundary['bins'], debit)
    delta = interval_ball(boundary['far_gram_trace_bound'], boundary['grid_bits'])
    if not sum(eta, arb(0)).overlaps(delta):
        raise ValueError('Changed Gamma trace-tail interface.')
    WU = arb_mat([[weights[i]*U[i, k] for k in range(U.ncols())] for i in range(M-N)])
    kap = interval_ball(old['kappa'], old['grid_bits'])
    last = selected[-1]
    data_s = {n: interval_ball(rows[n]['s'], scalar['grid_bits']) for n in range(N+1, M+1)}
    D0 = arb_mat(U.ncols(), len(selected))
    for first in range(N+1, 2*last+1, 128):
        modes = list(range(first, min(first+128, 2*last+1)))
        Z = arb_mat(len(modes), len(selected))
        for i, m in enumerate(modes):
            gm0, gm1 = arb((-1)**m), arb((-1)**m)*m*data_s[m]
            for j, n in enumerate(selected):
                Z[i, j] = (F[j, j] if m == n else
                    kap*(gm0*GE[j, 1]-gm1*GE[j, 0])/(m*m-n*n))
        D0 += arb_mat([[WU[m-N-1, k] for k in range(U.ncols())] for m in modes]).transpose()*Z
    print('new Gamma/support direct projection completed; seconds', time.monotonic()-started, flush=True)
    order = 384
    moments, absolute = arb_mat(order, 2*U.ncols()), arb_mat(1, 2*U.ncols())
    for first in range(2*last+1, M+1, 256):
        modes = list(range(first, min(first+256, M+1)))
        powers, values = arb_mat(len(modes), order), arb_mat(len(modes), 2*U.ncols())
        for i, m in enumerate(modes):
            ratio, power = arb(fmpq(last*last, m*m)), arb(1)
            for r in range(order):
                powers[i, r] = power
                power *= ratio
            gm = [arb((-1)**m), arb((-1)**m)*m*data_s[m]]
            for k in range(U.ncols()):
                for j in range(2):
                    value = WU[m-N-1, k]*gm[j]/(m*m)
                    values[i, 2*k+j] = value
                    absolute[0, 2*k+j] += abs(value)
        moments += powers.transpose()*values
    ratio = arb(fmpq(last, 2*last+1))
    remainder = ratio**(2*order)/(1-ratio*ratio)
    for k in range(U.ncols()):
        pa = arb_poly([moments[r, 2*k] for r in range(order)])
        pb = arb_poly([moments[r, 2*k+1] for r in range(order)])
        for j, n in enumerate(selected):
            t = arb(fmpq(n*n, last*last))
            error = kap*remainder*(abs(GE[j, 1])*absolute[0, 2*k]
                                   +abs(GE[j, 0])*absolute[0, 2*k+1])
            D0[k, j] += kap*(GE[j, 1]*pa(t)-GE[j, 0]*pb(t))+arb(0, error.upper())
    result = dict(schema='ACTUAL_C27_GAMMA_SUPPORT_ACTION_V1', metadata=provenance,
                  grid_bits=args.output_bits, working_bits=args.precision,
                  mixed_inverse_generator=matrix_out(inverse_mixed, args.output_bits),
                  finite_gamma_support=matrix_out(D0, args.output_bits),
                  direct_last_mode=2*last, geometric_order=order,
                  geometric_ratio=envelope(ratio, args.output_bits),
                  geometric_factor=envelope(remainder, args.output_bits),
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started)
    write_packet(args.out, result)
    print('new Gamma/support thin action archived; seconds', time.monotonic()-started, flush=True)


def consume(args):
    started = time.monotonic()
    old, metadata_bits, _, _, boundary, selected, GE, TZ, F, gram, inverse_GE, provenance = setup(args)
    saved, action_hash = read_json(args.finite)
    candidate, candidate_hash = read_json(args.direction)
    if saved['schema'] != 'ACTUAL_C27_GAMMA_SUPPORT_ACTION_V1' or saved['metadata'] != provenance:
        raise ValueError('Wrong new Gamma/support thin action.')
    if (candidate['schema'], candidate['original_target_sha256'], candidate['finite_action_sha256'],
            candidate['basis']) != ('ACTUAL_C27_JOINT_GAMMA_DIRECTION_V1', provenance['target_sha256'],
                                    provenance['boundary_sha256'], provenance['basis']):
        raise ValueError('Require the saved complete new negative direction.')
    vq = [fmpq(int(x), 1 << candidate['direction_dyadic_bits'])
          for x in candidate['direction_dyadic_integers']]
    if len(vq) != N+1:
        raise ValueError('Incomplete original retained vector.')
    v = arb_mat([[arb(x)] for x in vq])
    norm_v = arb(sum((x*x for x in vq), fmpq(0)))
    D0 = matrix_in(saved['finite_gamma_support'], saved['grid_bits'])
    inverse_mixed = matrix_in(saved['mixed_inverse_generator'], saved['grid_bits'])
    GP = arb_mat([[interval_ball(x, old['grid_bits']) for x in row[:2]] for row in old['generators']])
    TR = arb_mat([[interval_ball(x, old['grid_bits']) for x in row[2:4]] for row in old['generators']])
    kap = interval_ball(old['kappa'], old['grid_bits'])
    base_x, transformed_cross = arb_mat(768, 1), arb_mat(768, 1)
    for i, m in enumerate(selected):
        for n in range(N+1):
            base_x[i, 0] += kap*(inverse_GE[i, 0]*GP[n, 1]-inverse_GE[i, 1]*GP[n, 0])/(m*m-n*n)*v[n, 0]
            transformed_cross[i, 0] += kap*(inverse_mixed[i, 0]*GP[n, 1]-inverse_mixed[i, 1]*GP[n, 0]
                +inverse_GE[i, 0]*TR[n, 1]-inverse_GE[i, 1]*TR[n, 0])/(m*m-n*n)*v[n, 0]
    finite_gradient = F*transformed_cross
    Kfinite = gram-F
    scalar, rows, pilot, far = inputs()
    alpha, scale, Pq, Hq, _ = reconstruct_factors(dict(far, N=selected[-1]),
                                dict(scalar, rows=[rows[n] for n in range(selected[-1]+1)]))
    PE = arb_mat([[arb(Pq[i, n]) for n in selected] for i in range(Pq.nrows())])
    core = arb(alpha*scale)*arb_mat(Hq)
    PR = arb_mat([[interval_ball(row[4+k], old['grid_bits']) for row in old['generators']]
                 for k in range(core.nrows())])
    pr = PR*v
    epsilon = interval_ball(old['metadata']['far_matrix_error'], metadata_bits).upper()
    allowance = interval_ball(old['metadata']['complete_far_operator_allowance'], metadata_bits)
    qold = interval_ball(candidate['old_lower_target_quadratic'], candidate['grid_bits'])
    A0 = matrix_in(boundary['finite_gram'], boundary['grid_bits'])
    Beta0 = matrix_in(boundary['finite_beta'], boundary['grid_bits'])
    delta = interval_ball(boundary['far_gram_trace_bound'], boundary['grid_bits']).upper()
    refinement_path = getattr(args, 'channel_refinement', None)
    if refinement_path:
        refinement, refinement_hash = read_json(refinement_path)
        if (refinement['schema'], refinement['metadata'], refinement['parent_support_action_sha256']) != (
                'ACTUAL_C27_GAMMA_RATE_CONTRAST_V1', provenance, action_hash):
            raise ValueError('Rate contrasts must extend these exact parent channels and inputs.')
        rb = refinement['grid_bits']
        across = matrix_in(refinement['finite_cross'], rb)
        extra_A = matrix_in(refinement['finite_gram'], rb)
        extra_Beta = matrix_in(refinement['finite_beta'], rb)
        extra_D = matrix_in(refinement['finite_support'], rb)
        parent_count, extra_count = A0.nrows(), extra_A.nrows()
        if ((across.nrows(), across.ncols()), (extra_A.nrows(), extra_A.ncols()),
                (extra_Beta.nrows(), extra_Beta.ncols()), (extra_D.nrows(), extra_D.ncols())) != (
                (parent_count, extra_count), (extra_count, extra_count),
                (extra_count, N+1), (extra_count, len(selected))):
            raise ValueError('Incomplete common rate-contrast actions.')
        combined = arb_mat(parent_count+extra_count, parent_count+extra_count)
        for i in range(parent_count):
            for j in range(parent_count):
                combined[i, j] = A0[i, j]
            for j in range(extra_count):
                combined[i, parent_count+j] = across[i, j]
                combined[parent_count+j, i] = across[i, j]
        for i in range(extra_count):
            for j in range(extra_count):
                combined[parent_count+i, parent_count+j] = extra_A[i, j]
        A0 = combined
        Beta0 = arb_mat([[Beta0[i, j] if i < parent_count else extra_Beta[i-parent_count, j]
                         for j in range(N+1)] for i in range(parent_count+extra_count)])
        D0 = arb_mat([[D0[i, j] if i < parent_count else extra_D[i-parent_count, j]
                      for j in range(len(selected))] for i in range(parent_count+extra_count)])
        delta += interval_ball(refinement['far_gram_trace_bound'], rb).upper()
        provenance = dict(provenance, channel_refinement_sha256=refinement_hash)
    keep_even_pole = getattr(args, 'keep_even_pole', False)
    if keep_even_pole:
        if boundary['bins'][0] != [0, 1]:
            raise ValueError('The retained even pole must share the singleton b0=1/2 channel.')
        width = arb(27).log()
        pole_scale = (width/4).exp()
        for j in range(A0.ncols()):
            A0[0, j] *= pole_scale
        for i in range(A0.nrows()):
            A0[i, 0] *= pole_scale
        for j in range(Beta0.ncols()):
            Beta0[0, j] *= pole_scale
        for j in range(D0.ncols()):
            D0[0, j] *= pole_scale
        debit = arb(fmpq(int(pilot['far_weight_numerator']), int(pilot['far_weight_denominator'])))
        c0 = 1-(-width/2).exp()
        eta0 = debit*c0/width*(width/(2*arb.pi()))**4/(3*M**3)
        delta += (pole_scale**2-1)*eta0.upper()
        provenance = dict(provenance, retained_even_pole=True)
    beta_base = Beta0*v
    channel_count = A0.nrows()
    S0 = arb_mat(A0)
    for i in range(S0.nrows()):
        S0[i, i] += 1
    S0invD = S0.solve(D0, algorithm='precond')
    floor = fmpq(1, 16)
    S0inverse = S0.inv()
    floor_generators = arb_mat([[GE[i, j] for j in range(2)]+[TZ[i, j] for j in range(2)]
        +[D0[j, i] for j in range(channel_count)]+[m*m*D0[j, i] for j in range(channel_count)]
        for i, m in enumerate(selected)])
    floor_diagonal = [Kfinite[i, i]-sum((D0[j, i]*S0invD[j, i] for j in range(channel_count)), arb(0))-arb(floor)
                      for i in range(768)]
    floor_verdict, failure, floor_pivots, floor_order = eliminate(
        floor_diagonal, floor_generators, S0inverse, -kap, [m*m for m in selected],
        pivoting='max-diagonal')
    if floor_verdict != 'positive-pivots':
        raise ArithmeticError(f'Changed finite Gamma K floor not certified: {floor_verdict} at {failure}.')
    # Dropping the far squared residual before minimizing over channel z gives
    # K_G >= GramZ0-F-D0* S0^-1 D0, so these finite pivots bound every trial.
    trial_cache_hash = None
    if args.trial:
        previous, trial_cache_hash = read_json(args.trial)
        if (previous['schema'], previous['metadata'], previous['finite_action_sha256'],
                previous['direction_sha256']) != ('ACTUAL_C27_GAMMA_SUPPORT_OPTIMUM_V1',
                provenance, action_hash, candidate_hash):
            raise ValueError('Require a trial cache for these exact inputs and direction.')
        integers = [int(x) for x in previous['trial_correction_integers']]
        trial_bits = previous['trial_correction_bits']
        if len(integers) != 768 or trial_bits < 64:
            raise ValueError('Incomplete or invalid cached exact trial.')
    else:
        Kmid = Kfinite+PE.transpose()*core*PE-D0.transpose()*S0invD
        gmid = finite_gradient+PE.transpose()*core*pr-D0.transpose()*S0.solve(beta_base, algorithm='precond')
        proposed = Kmid.solve(gmid, algorithm='precond')
        trial_bits = args.trial_bits
        integers = [trial_integer(proposed[i, 0], trial_bits) for i in range(768)]
    t = arb_mat([[arb(fmpq(x, 1 << trial_bits))] for x in integers])
    u_norm_sq = norm_v+sum(((base_x[i, 0]+t[i, 0])**2 for i in range(768)), arb(0))
    pu = pr-PE*t
    far_center = (pu.transpose()*core*pu)[0, 0]
    far_error = epsilon*u_norm_sq.upper()
    old_center = (pr.transpose()*core*pr)[0, 0]
    finite_change = (2*t.transpose()*finite_gradient-t.transpose()*Kfinite*t)[0, 0]
    margin_center = qold+allowance*norm_v+finite_change+old_center-far_center
    margin_W = margin_center+arb(0, far_error.upper())
    far_upper = far_center+far_error
    if not far_upper > 0:
        raise ArithmeticError('Require a positive upper bound on actual far residual energy.')
    beta0 = beta_base-D0*t
    z0 = S0.solve(beta0, algorithm='precond')
    zintegers = [trial_integer(z0[i, 0], args.trial_bits) for i in range(S0.nrows())]
    z = arb_mat([[arb(fmpq(x, 1 << args.trial_bits))] for x in zintegers])
    norm = lambda x: sum((v*v for v in x.entries()), arb(0)).sqrt()
    e = (delta*far_upper.upper()).sqrt()
    gain_lower = (2*z.transpose()*beta0-z.transpose()*S0*z)[0, 0]-2*norm(z)*e-delta*norm(z)**2
    gain_upper = (beta0.transpose()*z0)[0, 0]+2*norm(beta0)*e+e**2
    lower, upper = margin_W+gain_lower, margin_W+gain_upper
    # The complete actual far support Gram, including its error, bounds D_infinity.
    far_support = PE.transpose()*core*PE
    trace_upper = sum((far_support[i, i] for i in range(768)), arb(0))+768*epsilon
    if not trace_upper > 0:
        raise ArithmeticError('Require a positive far support trace bound.')
    Dfar_norm = (delta*trace_upper.upper()).sqrt()
    z_error = e+delta*norm(z0)
    gradient = finite_gradient-Kfinite*t+PE.transpose()*core*pu-D0.transpose()*z0
    gradient_error = arb(768).sqrt()*epsilon*u_norm_sq.upper().sqrt()+norm(D0)*z_error+Dfar_norm*(norm(z0)+z_error)
    correction = (norm(gradient)+gradient_error)**2/arb(floor)
    optimal_upper = upper+correction
    verdict = ('all-trials-on-this-Gamma-support-excluded' if optimal_upper < 0 else
               'this-direction-repaired' if lower > 0 else 'inconclusive')
    result = dict(schema='ACTUAL_C27_GAMMA_SUPPORT_OPTIMUM_V1', metadata=provenance,
                  finite_action_sha256=action_hash, direction_sha256=candidate_hash,
                  grid_bits=args.output_bits, working_bits=args.precision,
                  trial_correction_bits=trial_bits, trial_correction_integers=list(map(str, integers)),
                  trial_cache_sha256=trial_cache_hash,
                  channel_count=channel_count,
                  retained_even_pole=keep_even_pole,
                  channel_trial_bits=args.trial_bits,
                  channel_trial_integers=list(map(str, zintegers)),
                  trial_base='Same exact Y0=E F^-1 B applied to the complete saved dyadic direction.',
                  trial_parameterization='y=Y0 v+E t; t is the saved exact dyadic correction.',
                  combined_norm_sq_upper=envelope(u_norm_sq.upper(), args.output_bits),
                  changed_K_floor=str(floor),
                  changed_K_floor_method='Drop the nonnegative far channel residual, then apply the existing structured LDL to GramZ0-F-D0* S0^-1 D0-I/16.',
                  changed_K_floor_positive_pivot_count=len(floor_pivots),
                  changed_K_floor_pivot_order=[selected[i] for i in floor_order],
                  changed_K_floor_pivot_intervals=[envelope(p, args.output_bits) for p in floor_pivots],
                  original_metadata_grid_bits=metadata_bits,
                  original_margin=envelope(qold, args.output_bits),
                  diagonal_margin=envelope(margin_W, args.output_bits),
                  far_residual_energy_upper=envelope(far_upper, args.output_bits),
                  Gamma_gain_lower=envelope(gain_lower, args.output_bits),
                  Gamma_gain_upper=envelope(gain_upper, args.output_bits),
                  updated_margin_lower=envelope(lower, args.output_bits), updated_margin_upper=envelope(upper, args.output_bits),
                  gradient_center_norm=envelope(norm(gradient), args.output_bits),
                  gradient_error_norm_upper=envelope(gradient_error, args.output_bits),
                  optimal_correction_upper=envelope(correction, args.output_bits),
                  all_trial_margin_upper=envelope(optimal_upper, args.output_bits),
                  verdict=verdict, elapsed_seconds=time.monotonic()-started,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  conditions=['Same actual common form/domain, exterior and original orthonormal mass.',
                    'Existing Gamma boundary/diagonal inverse comparison and sound saved full target and extended-support far enclosure.',
                    'Finite Gamma shorting bound (HE)*W_G(HE)-E*HE>=GramZ0-F-D0* S0^-1 D0>=I/16, with all768 directed pivots positive.'],
                  scope=f'Changed {channel_count}-channel W_G, same selected768-mode support and one complete saved retained direction. '
                    'A negative optimal upper excludes all trials on this support for this W_G only; a positive lower repairs this direction only. '
                    'Neither sign determines the actual Schur form, cofinal positivity or RH. No old scientific producer or old fixed-W optimizer rerun; no Lean verification.')
    args.out.write_text(json.dumps(result, indent=2)+'\n')
    print(verdict, 'margin bounds', lower.str(24), upper.str(24),
          'optimal upper', optimal_upper.str(24), 'seconds', time.monotonic()-started, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--produce', action='store_true')
    p.add_argument('--target', type=Path, default=HERE/'structured-joint-2048-target.json.zlib')
    p.add_argument('--pivots', type=Path, default=HERE/'structured-joint-2048-pivots.json')
    p.add_argument('--support', type=Path, default=HERE/'support-2048-finite-action.json.zlib')
    p.add_argument('--generators', type=Path, default=HERE/'joint-residual-generators-2048.json.zlib')
    p.add_argument('--boundary', type=Path, default=HERE/'joint-gamma-6765-finite.json.zlib')
    p.add_argument('--direction', type=Path, default=HERE/'joint-gamma-6765-direction.json')
    p.add_argument('--finite', type=Path)
    p.add_argument('--trial', type=Path, help='Reuse a matching saved exact trial without another proposal solve.')
    p.add_argument('--channel-refinement', type=Path, help='Append cached orthogonal rate contrasts to the same parent channels.')
    p.add_argument('--keep-even-pole', action='store_true', help='Retain the actual even pole by rescaling the singleton b0=1/2 channel.')
    p.add_argument('--precision', type=int, default=1024)
    p.add_argument('--output-bits', type=int, default=640)
    p.add_argument('--trial-bits', type=int, default=384)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if a.out.exists() or a.precision < max(640, a.output_bits, a.trial_bits)+192 or (not a.produce and a.finite is None):
        p.error('Use a fresh output,192 guard bits and a saved finite action unless producing it.')
    ctx.prec = a.precision
    produce(a) if a.produce else consume(a)


if __name__ == '__main__':
    main()
