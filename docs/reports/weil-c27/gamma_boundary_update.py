"""Retain finite positive even Gamma boundary channels in the inverse majorant.

Consume the saved full dyadic direction and support optimizer. Archive only
their previously unsaved residual action; then reuse it for thin boundary
projections. Existing arithmetic, diagonal weights and far energy are inputs.
No old target, optimizer, scalar or far-moment producer is run.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import zlib

import flint
from flint import acb, arb, arb_mat, arb_poly, ctx, fmpq

from arb_arithmetic import Arithmetic
from arb_inverse_residual_obstruction import HERE, M, N, inputs
from arb_inverse_residual_pilot import trial_integer
from diagonal_floor_witness import interval_ball
from structured_inverse_residual_interval import envelope


def read_json(path):
    raw = path.read_bytes()
    return json.loads(zlib.decompress(raw) if path.suffix == '.zlib' else raw), hashlib.sha256(raw).hexdigest()


def setup(args):
    candidate, ch = read_json(args.candidate)
    saved, sh = read_json(args.optimizer)
    if (candidate['schema'], saved['schema'], saved['candidate_sha256']) != (
            'ACTUAL_C27_SUPPORT_DIRECTION_V1', 'ACTUAL_C27_OPTIMAL_SUPPORT_TEST_V1', ch):
        raise ValueError('Require the saved matching complete direction and optimizer.')
    if (candidate['N'], candidate['M'], saved['N'], saved['M']) != (N, M, N, M):
        raise ValueError('Wrong common retained and middle ranges.')
    v = [fmpq(int(x), 1 << candidate['candidate_dyadic_bits'])
         for x in candidate['candidate_dyadic_integers']]
    x = [fmpq(int(x), 1 << saved['trial_dyadic_bits']) for x in saved['trial_dyadic_integers']]
    if len(v) != N+1 or len(x) != saved['last_mode']-N:
        raise ValueError('Incomplete source or trial.')
    scalar, rows, pilot, _ = inputs()
    arithmetic = Arithmetic(27, scalar['grid_bits'])
    weights = [arb(fmpq(int(x), 1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    D = fmpq(int(pilot['far_weight_numerator']), int(pilot['far_weight_denominator']))
    if len(weights) != M-N or not D > 0 or not all(w > 0 for w in weights):
        raise ValueError('Incomplete or nonpositive diagonal majorant.')
    provenance = dict(candidate_sha256=ch, optimizer_sha256=sh, N=N, M=M,
                      last_mode=saved['last_mode'], basis=saved['basis'],
                      input_sha256=saved['input_sha256'])
    return saved, v+[-t for t in x], scalar, rows, weights, arb(D), arithmetic, provenance


def prepare(args):
    started = time.monotonic()
    saved, uq, scalar, rows, weights, _, arithmetic, provenance = setup(args)
    last = len(uq)-1
    data = {n: (interval_ball(r['s'], scalar['grid_bits']), interval_ball(r['a'], scalar['grid_bits']))
            for n, r in rows.items()}
    u = [arb(x) for x in uq]
    g0 = [arb((-1)**n)/(arb(2).sqrt() if n == 0 else 1) for n in range(last+1)]
    g1 = [arb((-1)**n)*n*data[n][0] for n in range(last+1)]
    a, b = [g1[n]*u[n] for n in range(last+1)], [g0[n]*u[n] for n in range(last+1)]
    direct_end, order = 2*last, 128
    A, B = [], []
    power = [arb(1)]*(last+1)
    ratios = [arb(fmpq(n*n, last*last)) for n in range(last+1)]
    for _ in range(order):
        A.append(sum((a[n]*power[n] for n in range(last+1)), arb(0)))
        B.append(sum((b[n]*power[n] for n in range(last+1)), arb(0)))
        power = [power[n]*ratios[n] for n in range(last+1)]
    pa, pb = arb_poly(A), arb_poly(B)
    sum_a, sum_b = sum(map(abs, a), arb(0)), sum(map(abs, b), arb(0))
    kappa = 2/arithmetic.pi
    residual = []
    sample_checks = []
    for m in range(N+1, M+1):
        if m <= direct_end:
            value = sum((arithmetic.even_entry(m,n,data)*u[n] for n in range(last+1) if uq[n]), arb(0))
        else:
            t = arb(fmpq(last*last, m*m))
            value = kappa*((-1)**m)/(m*m)*(pa(t)-m*data[m][0]*pb(t))
            radius = kappa/(m*m)*t**order/(1-t)*(sum_a+m*abs(data[m][0])*sum_b)
            value += arb(0, radius.upper())
            if m in (direct_end+1, 8192, 16384, M):
                direct = sum((arithmetic.even_entry(m,n,data)*u[n] for n in range(last+1) if uq[n]), arb(0))
                if not value.overlaps(direct):
                    raise ArithmeticError('Residual geometric expansion disagrees with actual direct row.')
                sample_checks.append(m)
        residual.append(value)
        if (m-N) % 4096 == 0:
            print('new residual action through',m,'seconds',time.monotonic()-started,flush=True)
    energy = sum((w*r*r for w,r in zip(weights,residual)), arb(0))
    old_energy = interval_ball(saved['middle_residual_energy'], saved['grid_bits'])
    if not energy.overlaps(old_energy):
        raise ArithmeticError('Residual action does not reproduce the archived finite energy.')
    result = dict(schema='ACTUAL_C27_RESIDUAL_ACTION_V1', metadata=provenance,
                  grid_bits=args.output_bits, working_bits=args.precision,
                  residual=[envelope(r,args.output_bits) for r in residual],
                  middle_energy=envelope(energy,args.output_bits),
                  direct_last_mode=direct_end, geometric_order=order,
                  direct_geometric_sample_checks=sample_checks,
                  producer_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started)
    args.out.write_bytes(zlib.compress(json.dumps(result,separators=(',',':')).encode(),9))
    print('residual archived; finite energy',energy.str(20),'seconds',time.monotonic()-started,flush=True)


def fibonacci_bins(last):
    if last < 8:
        raise ValueError('Use last-bin>=8.')
    endpoints = [0,1,2,3,5,8]
    while endpoints[-1] < last:
        endpoints.append(min(last,endpoints[-1]+endpoints[-2]))
    return list(zip(endpoints,endpoints[1:]))


def boundary_columns(arithmetic, bins, D):
    """Actual finite boundary columns and their existing trace-tail allowance."""
    c = [1-(-arb(fmpq(4*a+1,2))*arithmetic.d).exp() for a,_ in bins]
    factors = [(c[k]/(b-a)/arithmetic.d).sqrt() for k,(a,b) in enumerate(bins)]
    U = arb_mat(M-N,len(bins))
    for i,m in enumerate(range(N+1,M+1)):
        im = arithmetic.pi*m/arithmetic.d
        values = {e:acb(arb(e)+arb(fmpq(1,4)),im).digamma().real
                  for e in {x for ab in bins for x in ab}}
        for k,(a,b) in enumerate(bins):
            U[i,k] = (-1)**m*factors[k]*(values[b]-values[a])
    # Each column has |U_m|<=a_k/m^2. The far Gram is bounded by its trace,
    # not by an entrywise outer product; mixed tails use weighted Cauchy.
    constants = [2*factors[k]*(b-a)*arb(fmpq(2*a+2*b-1,2))
                 *(arithmetic.d/(2*arithmetic.pi))**2 for k,(a,b) in enumerate(bins)]
    tail_energies = [D*t*t/(3*M**3) for t in constants]
    return U,tail_energies


def certify(args):
    started = time.monotonic()
    saved, _, _, _, weights, D, arithmetic, provenance = setup(args)
    packet, rh = read_json(args.residual)
    if packet['schema'] != 'ACTUAL_C27_RESIDUAL_ACTION_V1' or packet['metadata'] != provenance:
        raise ValueError('Residual action belongs to another exact source/trial.')
    r = [interval_ball(x,packet['grid_bits']) for x in packet['residual']]
    if len(r) != M-N:
        raise ValueError('Require the complete actual finite residual.')
    bins = fibonacci_bins(args.last_bin)
    U,tail_energies = boundary_columns(arithmetic,bins,D)
    print('new thin boundary projections assembled; seconds',time.monotonic()-started,flush=True)
    WU = arb_mat([[weights[i]*U[i,k] for k in range(len(bins))] for i in range(M-N)])
    A = U.transpose()*WU
    beta = WU.transpose()*arb_mat([[x] for x in r])
    delta = sum(tail_energies,arb(0))
    far = interval_ball(saved['far_residual_energy'],saved['grid_bits'])
    if not far > 0:
        raise ArithmeticError('Require a positive archived actual residual far-energy enclosure.')
    beta_errors = [(t*far.upper()).sqrt() for t in tail_energies]
    S = A+arb_mat([[arb(1) if j == k else arb(0) for k in range(len(bins))]
                   for j in range(len(bins))])
    z0 = S.solve(beta,algorithm='precond')
    zintegers = [trial_integer(z0[k,0],args.trial_bits) for k in range(len(bins))]
    zq = [fmpq(x,1 << args.trial_bits) for x in zintegers]
    z = arb_mat([[arb(x)] for x in zq])
    midpoint_gain = (2*z.transpose()*beta-z.transpose()*S*z)[0,0]
    loss = 2*sum((abs(arb(zq[k]))*beta_errors[k] for k in range(len(bins))),arb(0))
    loss += delta*arb(sum((x*x for x in zq),fmpq(0)))
    gain_lower = midpoint_gain-loss
    margin = interval_ball(saved['fixed_trial_margin'],saved['grid_bits'])
    new_lower = margin+gain_lower
    # S_actual>=I+A_middle, so completing a second finite square gives an
    # upper gain including every admissible far mixed projection.
    beta_norm = sum((beta[k,0]**2 for k in range(len(bins))),arb(0)).sqrt()
    error_norm = sum((e*e for e in beta_errors),arb(0)).sqrt()
    gain_upper = (beta.transpose()*S.solve(beta,algorithm='precond'))[0,0]+2*beta_norm*error_norm+error_norm**2
    new_upper = margin+gain_upper
    verdict = ('saved-direction-repaired' if new_lower > 0 else
               'saved-trial-still-excluded' if new_upper < 0 else 'inconclusive')
    result = dict(schema='ACTUAL_C27_GAMMA_BOUNDARY_UPDATE_V1',metadata=provenance,
                  residual_sha256=rh,working_bits=args.precision,grid_bits=args.output_bits,
                  bins=bins,trial_dyadic_bits=args.trial_bits,trial_dyadic_integers=[str(x) for x in zintegers],
                  trial_rationals=[str(x) for x in zq],
                  tail_gram_trace_bound=envelope(delta,args.output_bits),
                  mixed_tail_norm_bound=envelope(error_norm,args.output_bits),
                  gain_lower=envelope(gain_lower,args.output_bits),gain_upper=envelope(gain_upper,args.output_bits),
                  original_margin=saved['fixed_trial_margin'],original_margin_grid_bits=saved['grid_bits'],
                  updated_margin_lower=envelope(new_lower,args.output_bits),updated_margin_upper=envelope(new_upper,args.output_bits),
                  verdict=verdict,elapsed_seconds=time.monotonic()-started,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  python_flint=flint.__version__,FLINT=flint.__FLINT_VERSION__,
                  conditions=['Same actual complete even common-form and exterior-domain identifications.',
                              'Archived diagonal W dominates the reciprocal even Gamma-minus-prime lower comparison.',
                              'Existing even Gamma boundary positive series; finite-bin Cauchy lower update.',
                              'Sound saved exact trial, original margin, actual finite residual and far energy.'],
                  scope='Retains actual infinite-support positive Gamma boundary vectors in the existing Woodbury inverse majorant. '
                        'Tests only the archived exact retained direction AND archived fixed support trial, not all trials, full Schur positivity or RH.')
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(verdict,'gain',gain_lower.str(20),gain_upper.str(20),'margin',new_lower.str(20),new_upper.str(20),
          'seconds',time.monotonic()-started,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepare',action='store_true')
    p.add_argument('--candidate',type=Path,default=HERE/'support-2048-direction.json')
    p.add_argument('--optimizer',type=Path,default=HERE/'support-2048-optimized.json')
    p.add_argument('--residual',type=Path)
    p.add_argument('--last-bin',type=int,default=6765)
    p.add_argument('--precision',type=int,default=512)
    p.add_argument('--output-bits',type=int,default=256)
    p.add_argument('--trial-bits',type=int,default=256)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists() or a.precision < max(a.output_bits,a.trial_bits)+128 or (not a.prepare and a.residual is None):
        p.error('Use a fresh output, at least128 guard bits and a saved residual unless preparing it.')
    ctx.prec=a.precision
    prepare(a) if a.prepare else certify(a)


if __name__ == '__main__':
    main()
