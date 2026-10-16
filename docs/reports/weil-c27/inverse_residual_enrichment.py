"""Propose and enclose actual-exterior inverse corrections on one fixed direction.

Reuse the published complete scalar packets, diagonal weights and far moments.
The retained direction is unchanged. A finite H solve proposes a correction;
only a separately enclosed full residual can accept an infinite-exterior bound.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time
from types import SimpleNamespace

import flint
from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

from arb_arithmetic import Arithmetic, grid_interval
from arb_far_coupling import reconstruct_factors, far_bounds
from arb_inverse_residual_obstruction import HERE, N, M, inputs, direction
from diagonal_floor_witness import interval_ball


BASELINE_SHA = '3b850c5e42f0c5671884a2abdabbc66b5e7611df9bfb6e24f63b512215679a51'


def baseline():
    path = HERE/'inverse-residual-k128-obstruction.json'
    if hashlib.sha256(path.read_bytes()).hexdigest() != BASELINE_SHA:
        raise ValueError('Use the published full-direction baseline.')
    return json.loads(path.read_text())


def propose(args):
    import numpy as np
    from threadpoolctl import threadpool_limits
    started = time.monotonic()
    scalar, rows, pilot, far = inputs()
    old = direction(args.candidate)
    base = baseline()
    if (old['columns'], hashlib.sha256(args.candidate.read_bytes()).hexdigest()) != (128, base['candidate_sha256']):
        raise ValueError('Use the same exact128-column retained direction/trial.')
    last = args.last_mode
    if not N+128 < last < M:
        raise ValueError('Enrichment must leave the excluded span and stay below M.')
    bits = scalar['grid_bits']
    def midpoint(row, key):
        return math.ldexp(int(row[key][0])+int(row[key][1]), -bits-1)
    s = np.array([midpoint(rows[n], 's') for n in range(M+1)])
    a = np.array([midpoint(rows[n], 'a') for n in range(last+1)])
    indices = np.arange(last+1)
    phases = np.where(indices % 2, -1., 1.)
    norms = np.ones(last+1)
    norms[0] = 1/math.sqrt(2)
    symbols = indices*s[:last+1]
    diagonal = a.copy()
    diagonal[1:] -= s[1:last+1]/(math.pi*indices[1:])
    def block(modes):
        den = modes[:, None]**2-indices[None, :]**2
        out = np.divide(symbols[None, :]-modes[:, None]*s[modes, None], den,
                        out=np.zeros_like(den, dtype=float), where=den != 0)
        rn = np.where(modes == 0, 1/math.sqrt(2), 1.)
        out *= (2/math.pi*np.where(modes % 2, -1., 1.)*rn)[:, None]*(phases*norms)[None, :]
        for i, m in enumerate(modes):
            if m <= last:
                out[i, m] = diagonal[m]
        return out
    v = np.ldexp(np.array([float(int(x)) for x in old['candidate_dyadic_integers']]), -old['candidate_dyadic_bits'])
    old_y = np.zeros(last-N)
    old_y[:128] = np.ldexp(np.array([float(int(x)) for x in old['trial_dyadic_integers']]), -old['trial_dyadic_bits'])
    weights = np.ldexp(np.array([float(int(x)) for x in pilot['middle_weight_numerators']]), -pilot['middle_weight_bits'])
    def old_value(key):
        pair = base[key]
        return math.ldexp(int(pair[0])+int(pair[1]), -base['grid_bits']-1)
    with threadpool_limits(limits=1):
        selected = np.arange(N+1, last+1)
        selected_rows = block(selected)
        H = selected_rows[:, N+1:]
        g = selected_rows[:, :N+1]@v
        r = g-H@old_y
        delta = np.linalg.solve(H, r)
        new_y = old_y+delta
        u = np.concatenate((v, -new_y))
        q_change = -2*delta@r+delta@H@delta
        q_new = old_value('finite_combined_Q')+q_change
        middle = 0.
        for first in range(N+1, M+1, 1024):
            modes = np.arange(first, min(first+1024, M+1))
            residual = block(modes)@u
            middle += residual@(weights[modes-N-1]*residual)
        alpha, D, P0, H0, _ = reconstruct_factors(
            dict(far, N=last), dict(scalar, rows=[rows[n] for n in range(last+1)]))
        P = np.array([[float(P0[i,j]) for j in range(last+1)] for i in range(P0.nrows())])
        Hm = np.array([[float(H0[i,j]) for j in range(H0.ncols())] for i in range(H0.nrows())])
        feature = P@u
        far_debit = float(alpha*D)*(feature@Hm@feature)
        result = {
            'schema':'ACTUAL_C27_INVERSE_ENRICHMENT_CANDIDATE_V1',
            'N':N, 'M':M, 'last_mode':last,
            'retained_candidate_sha256':hashlib.sha256(args.candidate.read_bytes()).hexdigest(),
            'baseline_sha256':BASELINE_SHA,
            'trial_dyadic_bits':96,
            'trial_dyadic_integers':[str(round(float(x)*(1 << 96))) for x in new_y],
            'proposal':{
                'algorithm':'Finite actual H solve for correction to old actual residual, not the old K optimizer.',
                'selected_residual_norm':float(np.linalg.norm(g-H@new_y)),
                'middle_residual_energy':float(middle),
                'far_residual_center':float(far_debit),
                'total_residual_energy':float(middle+far_debit),
                'combined_Q_from_baseline_delta':float(q_new),
                'margin_from_baseline_delta':float(q_new-middle-far_debit),
                'old_residual_energy':old_value('middle_residual_debit')+old_value('far_residual_debit'),
                'old_margin':old_value('fixed_trial_margin'),
                'elapsed_seconds':time.monotonic()-started,
            },
            'scope':'Floating proposal only; no directed sign or full-exterior inverse certificate. Same complete retained dyadic direction; original orthonormal mass.',
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        }
        args.out.write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(result['proposal'], indent=2), flush=True)


def certify(args):
    if args.precision < max(640, args.output_bits)+96 or args.output_bits < 64:
        raise ValueError('Keep at least96 working guard bits above the640-bit input.')
    started = time.monotonic()
    scalar, rows, pilot, far = inputs()
    old = direction(args.candidate)
    base = baseline()
    candidate = json.loads(args.trial.read_text())
    old_hash = hashlib.sha256(args.candidate.read_bytes()).hexdigest()
    if (candidate['schema'], candidate['N'], candidate['M'],
            candidate['retained_candidate_sha256'], candidate['baseline_sha256'],
            old['columns'], base['candidate_sha256']) != (
            'ACTUAL_C27_INVERSE_ENRICHMENT_CANDIDATE_V1', N, M,
            old_hash, BASELINE_SHA, 128, old_hash):
        raise ValueError('Require the unchanged complete direction and pinned baseline.')
    last = candidate['last_mode']
    if type(last) is not int or not N+128 < last < M:
        raise ValueError('Require an actual enrichment within scalar coverage.')
    yq = [fmpq(int(x), 1 << candidate['trial_dyadic_bits'])
          for x in candidate['trial_dyadic_integers']]
    if len(yq) != last-N:
        raise ValueError('Incomplete finite-domain trial.')
    vq = [fmpq(int(x), 1 << old['candidate_dyadic_bits'])
          for x in old['candidate_dyadic_integers']]
    old_yq = [fmpq(int(x), 1 << old['trial_dyadic_bits'])
              for x in old['trial_dyadic_integers']]+[fmpq(0)]*(last-N-128)
    dq = [a-b for a,b in zip(yq, old_yq)]
    uq = vq+[-x for x in yq]
    old_u = arb_mat([[arb(x)] for x in vq+[-x for x in old_yq]])
    delta = arb_mat([[arb(x)] for x in dq])
    u = arb_mat([[arb(x)] for x in uq])
    norm_u_sq = sum((x*x for x in uq), fmpq(0))
    norm_v_sq = sum((x*x for x in vq), fmpq(0))
    if not norm_v_sq > 0:
        raise ValueError('Nonzero complete retained direction required.')
    data = {n:(interval_ball(row['s'], scalar['grid_bits']),
               interval_ball(row['a'], scalar['grid_bits'])) for n,row in rows.items()}
    arithmetic = Arithmetic(27, scalar['grid_bits'])
    weights = [arb(fmpq(int(x), 1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    if len(weights) != M-N or not all(w > 0 for w in weights):
        raise ValueError('Require the complete positive archived weights.')

    # Reuse the old directed Q, evaluating only its exact correction.
    # The same selected rows also supply the first part of the new residual.
    selected_rows = arb_mat([[arithmetic.even_entry(m,n,data)
                             for n in range(last+1)] for m in range(N+1,last+1)])
    r_old = selected_rows*old_u
    Hdelta = arb_mat([[sum((selected_rows[i,N+1+j]*delta[j,0]
                            for j in range(last-N)), arb(0))]
                      for i in range(last-N)])
    q_change = -2*(delta.transpose()*r_old)[0,0]+(delta.transpose()*Hdelta)[0,0]
    q_new = interval_ball(base['finite_combined_Q'], base['grid_bits'])+q_change
    selected_residual = r_old-Hdelta
    selected_energy = sum((weights[i]*selected_residual[i,0]**2
                           for i in range(last-N)), arb(0))
    middle = selected_energy
    print('selected rows and Q correction enclosed; seconds', time.monotonic()-started, flush=True)
    for first in range(last+1, M+1, 128):
        modes = range(first, min(first+128, M+1))
        block = arb_mat([[arithmetic.even_entry(m,n,data)
                          for n in range(last+1)] for m in modes])
        residual = block*u
        middle += sum((weights[m-N-1]*residual[i,0]**2 for i,m in enumerate(modes)), arb(0))
        if (first-last-1) % 4096 == 0:
            print('middle through', modes.stop-1, 'seconds', time.monotonic()-started, flush=True)

    # The same168 moments apply, with new amplification through last.
    alpha, D, P, Hm, _ = reconstruct_factors(
        dict(far, N=last), dict(scalar, rows=[rows[n] for n in range(last+1)]))
    feature = P*fmpq_mat([[x] for x in uq])
    far_center = alpha*D*(feature.transpose()*Hm*feature)[0,0]
    symbol_bound = 2*arithmetic.H+sum((w for w,_ in arithmetic.primes), arb(0))
    scalar_radius = max(arb(int(rows[n]['s'][1])-int(rows[n]['s'][0]))
                        /(1 << (scalar['grid_bits']+1)) for n in range(last+1))
    bounds = far_bounds(SimpleNamespace(M=M, L=arithmetic.d, pi=arithmetic.pi,
                                       budget=symbol_bound),
                        last, far['jet_order'], scalar_radius,
                        far['grid_bits'], far['moment_diagonals'])
    epsilon = arb(D)*bounds['unweighted_matrix_error']
    far_debit = arb(far_center)+arb(0, (epsilon*arb(norm_u_sq)).upper())
    residual_energy = middle+far_debit
    lower_target = q_new-residual_energy
    old_residual = (interval_ball(base['middle_residual_debit'],base['grid_bits'])
                    +interval_ball(base['far_residual_debit'],base['grid_bits']))
    old_margin = interval_ball(base['fixed_trial_margin'],base['grid_bits'])
    improvement = lower_target-old_margin
    if not lower_target > 0 or not residual_energy < old_residual:
        raise ArithmeticError('No positive lower target with a reduced complete actual residual.')
    result = {
        'schema':'ACTUAL_C27_INVERSE_ENRICHMENT_RESULT_V1',
        'N':N, 'M':M, 'last_mode':last, 'basis':scalar['basis'],
        'working_bits':args.precision, 'grid_bits':args.output_bits,
        'candidate_sha256':hashlib.sha256(args.trial.read_bytes()).hexdigest(),
        'retained_candidate_sha256':old_hash, 'baseline_sha256':BASELINE_SHA,
        'retained_norm_sq':str(norm_v_sq), 'combined_norm_sq':str(norm_u_sq),
        'Q_correction':grid_interval(q_change,args.output_bits),
        'combined_Q':grid_interval(q_new,args.output_bits),
        'selected_residual_energy':grid_interval(selected_energy,args.output_bits),
        'middle_residual_energy':grid_interval(middle,args.output_bits),
        'far_residual_energy':grid_interval(far_debit,args.output_bits),
        'complete_residual_energy':grid_interval(residual_energy,args.output_bits),
        'lower_target':grid_interval(lower_target,args.output_bits),
        'improvement':grid_interval(improvement,args.output_bits),
        'extended_far_error_grid_bits':640,
        'extended_far_matrix_error':grid_interval(epsilon,640),
        'positive_directional_lower_target':True,
        'residual_strictly_reduced':True,
        'middle_rows':{'first':N+1,'last':M,'count':M-N,'selected_reused_in_sum':last-N},
        'elapsed_seconds':time.monotonic()-started,
        'python_flint':flint.__version__, 'FLINT':flint.__FLINT_VERSION__,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'dependency_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['arb_arithmetic.py','arb_far_coupling.py','arb_tail_moments.py',
             'arb_middle_scalars.py','arb_inverse_residual_pilot.py',
             'arb_inverse_residual_obstruction.py','diagonal_floor_witness.py']},
        'input_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['arithmetic-scalars.json.zlib','diagonal-floor-witness.json','middle-scalars.json.zlib',
             'far-coupling.json.zlib','inverse-residual-pilot.json.zlib']},
        'conditions':['Existing actual complete Fourier/form identification and exterior-domain premises.',
                      'The same positive actual exterior H has inverse majorized by the archived W.'],
        'scope':'Directed conditional scalar lower bound for the actual Schur quadratic on the unchanged full dyadic v128, by the existing inverse-residual identity. '
                'Finite H solving only proposes y. All middle/far residuals are paid; no finite compression is substituted for the infinite inverse. '
                'No full retained-matrix positivity, cofinal positivity, new general theorem, Lean verification or RH proof.',
    }
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print('complete residual',residual_energy.str(30),'lower target',lower_target.str(30),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path, default=HERE/'inverse-residual-k128-candidate.json')
    parser.add_argument('--last-mode', type=int, default=2048)
    parser.add_argument('--trial', type=Path, help='Enclose this fixed dyadic trial instead of proposing one.')
    parser.add_argument('--precision', type=int, default=768)
    parser.add_argument('--output-bits', type=int, default=384)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output path.')
    ctx.prec = args.precision
    if args.trial:
        certify(args)
    else:
        propose(args)


if __name__ == '__main__':
    main()
