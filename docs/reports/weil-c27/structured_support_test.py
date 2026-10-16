"""Test all trials on a fixed support along one saved-target Schur direction.

Reuse the complete directed target to select a finite dyadic direction, then
apply the existing square-completion upper bound. The finite weighted Gram
uses Trummer displacement and separate diagonals. No old scalar, moment,
full retained Gram or target producer runs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import zlib

import flint
from flint import arb, arb_mat, ctx, fmpq, fmpq_mat

from arb_arithmetic import Arithmetic
from arb_far_coupling import reconstruct_factors
from arb_inverse_residual_obstruction import HERE, M, N, inputs
from arb_inverse_residual_pilot import trial_integer
from diagonal_floor_witness import interval_ball
from structured_inverse_residual_interval import eliminate, envelope
from structured_residual_schur_finish import load_target, reconstruct_target


def prepare(args):
    started = time.monotonic()
    packet, prefix, remaining, target_hash, pivot_raw = load_target(args.target, args.pivots)
    mode = args.mode
    if mode not in remaining:
        raise ValueError('Use a coordinate in the saved unresolved Schur block.')
    target = reconstruct_target(packet)
    positive = arb_mat([[target[i,j] for j in prefix] for i in prefix])
    coupling = arb_mat([[target[i,mode]] for i in prefix])
    lift = positive.solve(coupling, algorithm='precond')
    integers = [0]*(N+1)
    integers[mode] = 1 << args.direction_bits
    for i,n in enumerate(prefix):
        integers[n] = trial_integer(-lift[i,0],args.direction_bits)
    exact = [fmpq(x,1 << args.direction_bits) for x in integers]
    v = arb_mat([[arb(x)] for x in exact])
    quadratic = (v.transpose()*target*v)[0,0]
    if not quadratic < 0:
        raise ArithmeticError('Rounded lifted direction did not retain a negative target quadratic.')
    result = {
        'schema':'ACTUAL_C27_SUPPORT_DIRECTION_V1','N':N,'M':M,
        'last_mode':packet['metadata']['last_mode'],'basis':packet['metadata']['basis'],
        'candidate_dyadic_bits':args.direction_bits,
        'candidate_dyadic_integers':[str(x) for x in integers],
        'candidate_norm_sq':str(sum((x*x for x in exact),fmpq(0))),
        'target_packet_sha256':target_hash,
        'pivot_result_sha256':hashlib.sha256(pivot_raw).hexdigest(),
        'selected_remaining_mode':mode,'grid_bits':args.output_bits,
        'target_quadratic':envelope(quadratic,args.output_bits),
        'working_bits':args.precision,'elapsed_seconds':time.monotonic()-started,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Exact finite dyadic retained direction from the saved lower-target Schur lift. '
                'Negative target energy is not negative actual Weil energy or exclusion of all support trials.',
    }
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print('dyadic target quadratic',quadratic.str(24),'seconds',time.monotonic()-started,flush=True)


def gram_from_generators(G,T,diagonal,kappa,modes):
    """Same finite weighted-Gram product rule; diagonal is extra input."""
    size = len(modes)
    gram = arb_mat(size,size)
    for i,m in enumerate(modes):
        gram[i,i] = diagonal[i]
        for j,n in enumerate(modes[:i]):
            numerator = (T[i,0]*G[j,1]-T[i,1]*G[j,0]
                         +G[i,0]*T[j,1]-G[i,1]*T[j,0])
            value = kappa*numerator/(m*m-n*n)
            gram[i,j] = value
            gram[j,i] = value
    return gram


def self_check():
    # Overlapping row/column sets require the actual diagonal, unlike a
    # disjoint cross block. Compare against an independent exact dense Gram.
    selected, outside = [2,4,7], [2,4,7,9,11]
    kap = fmpq(2,3)
    def gen(n):
        return [fmpq((-1)**n),fmpq((-1)**n*(n+1),n+2)]
    G = fmpq_mat([gen(n) for n in selected])
    GO = fmpq_mat([gen(n) for n in outside])
    Z = fmpq_mat([[fmpq(50+n) if m == n else
          kap*(gen(m)[0]*gen(n)[1]-gen(m)[1]*gen(n)[0])/(m*m-n*n)
          for n in selected] for m in outside])
    weights = [fmpq(i+2,i+1) for i in range(len(outside))]
    W = fmpq_mat([[weights[i] if i == j else 0 for j in range(len(outside))]
                  for i in range(len(outside))])
    T = Z.transpose()*W*GO
    exact = Z.transpose()*W*Z
    reconstructed = gram_from_generators(arb_mat(G),arb_mat(T),
                    [arb(exact[i,i]) for i in range(len(selected))],arb(kap),selected)
    assert all(reconstructed[i,j].contains(arb(exact[i,j]))
               for i in range(len(selected)) for j in range(len(selected)))
    F = fmpq_mat([[Z[i,j] for j in range(len(selected))] for i in range(len(selected))])
    K = exact-F
    generators = arb_mat([[arb(G[i,0]),arb(G[i,1]),arb(T[i,0]),arb(T[i,1])]
                           for i in range(len(selected))])
    verdict,_,pivots,order = eliminate([arb(K[i,i]-1) for i in range(len(selected))],
                generators,arb_mat(0,0),arb(-kap),[n*n for n in selected])
    assert verdict == 'positive-pivots'
    dense = fmpq_mat([[K[i,j]-(1 if i == j else 0) for j in range(len(selected))]
                      for i in range(len(selected))])
    for i in range(len(selected)):
        assert pivots[i].contains(arb(dense[i,i]))
        for j in range(i+1,len(selected)):
            for k in range(i+1,len(selected)):
                dense[j,k] -= dense[j,i]*dense[i,k]/dense[i,i]
    print('overlapping-row exact weighted Gram and signed K-floor LDL checks passed',flush=True)


def finite_action(args, selected, vq, scalar, data, weights, arithmetic, provenance):
    """Archive expensive finite sums once; later floor/optimizer consumers reuse them."""
    started = time.monotonic()
    count = len(selected)
    if args.finite:
        raw = args.finite.read_bytes()
        saved = json.loads(zlib.decompress(raw))
        if saved['schema'] != 'ACTUAL_C27_SUPPORT_ACTION_V1' or saved['metadata'] != provenance:
            raise ValueError('Finite support action belongs to another direction/target/input.')
    else:
        v = arb_mat([[arb(x)] for x in vq])
        def generator(n):
            phase = arb(-1 if n % 2 else 1)/(arb(2).sqrt() if n == 0 else 1)
            return [phase,phase*n*data[n][0]]
        GE = arb_mat([generator(n) for n in selected])
        Fdiagonal = [arithmetic.even_entry(n,n,data) for n in selected]
        Bv = arb_mat([[arithmetic.even_entry(m,n,data) for n in range(N+1)] for m in selected])*v
        Av = arb_mat([[arithmetic.even_entry(m,n,data) for n in range(N+1)] for m in range(N+1)])*v
        qv = (v.transpose()*Av)[0,0]
        print('new direction and selected finite block assembled; seconds',time.monotonic()-started,flush=True)
        T = arb_mat(count,2)
        diagonal = [arb(0)]*count
        weighted_cv = arb_mat(count,1)
        cv_energy = arb(0)
        for first in range(N+1,M+1,128):
            modes = list(range(first,min(first+128,M+1)))
            Z = arb_mat([[arithmetic.even_entry(m,n,data) for n in selected] for m in modes])
            Cv = arb_mat([[arithmetic.even_entry(m,n,data) for n in range(N+1)] for m in modes])*v
            WG = arb_mat([[weights[m-N-1]*g for g in generator(m)] for m in modes])
            WCv = arb_mat([[weights[m-N-1]*Cv[i,0]] for i,m in enumerate(modes)])
            T += Z.transpose()*WG
            weighted_cv += Z.transpose()*WCv
            cv_energy += (Cv.transpose()*WCv)[0,0]
            for j in range(count):
                diagonal[j] += sum((weights[m-N-1]*Z[i,j]**2 for i,m in enumerate(modes)),arb(0))
            if (first-N-1) % 4096 == 0:
                print('new support Gram/action through',modes[-1],'seconds',time.monotonic()-started,flush=True)
        bits = 640
        def matrix_out(A):
            return [[envelope(A[i,j],bits) for j in range(A.ncols())] for i in range(A.nrows())]
        saved = {
            'schema':'ACTUAL_C27_SUPPORT_ACTION_V1','metadata':provenance,
            'grid_bits':bits,'producer_working_bits':args.precision,
            'GE':matrix_out(GE),'T':matrix_out(T),
            'F_diagonal':[envelope(x,bits) for x in Fdiagonal],
            'gram_diagonal':[envelope(x,bits) for x in diagonal],
            'Bv':matrix_out(Bv),'weighted_cv':matrix_out(weighted_cv),
            'qv':envelope(qv,bits),'cv_energy':envelope(cv_energy,bits),
        }
        raw = zlib.compress(json.dumps(saved,separators=(',',':')).encode(),9)
        if args.finite_out:
            args.finite_out.write_bytes(raw)
            print('finite support action archived; bytes',len(raw),flush=True)
    bits = saved['grid_bits']
    def matrix_in(rows):
        return arb_mat([[interval_ball(x,bits) for x in row] for row in rows])
    GE,T = matrix_in(saved['GE']),matrix_in(saved['T'])
    if (GE.nrows(),GE.ncols(),T.nrows(),T.ncols()) != (count,2,count,2):
        raise ValueError('Incomplete finite generators.')
    Bv,weighted_cv = matrix_in(saved['Bv']),matrix_in(saved['weighted_cv'])
    if (Bv.nrows(),Bv.ncols(),weighted_cv.nrows(),weighted_cv.ncols()) != (count,1,count,1):
        raise ValueError('Incomplete finite actions.')
    if len(saved['F_diagonal']) != count or len(saved['gram_diagonal']) != count:
        raise ValueError('Incomplete independent finite diagonals.')
    F = arb_mat(count,count)
    kappa = 2/arithmetic.pi
    for i,m in enumerate(selected):
        F[i,i] = interval_ball(saved['F_diagonal'][i],bits)
        for j,n in enumerate(selected[:i]):
            F[i,j] = kappa*(GE[i,0]*GE[j,1]-GE[i,1]*GE[j,0])/(m*m-n*n)
            F[j,i] = F[i,j]
    diagonal = [interval_ball(x,bits) for x in saved['gram_diagonal']]
    return GE,T,F,Bv,weighted_cv,diagonal,interval_ball(saved['qv'],bits),interval_ball(saved['cv_energy'],bits),hashlib.sha256(raw).hexdigest()


def certify(args):
    started = time.monotonic()
    packet,_,_,target_hash,pivot_raw = load_target(args.target,args.pivots)
    old_pivots = json.loads(pivot_raw)
    candidate_raw = args.candidate.read_bytes()
    candidate = json.loads(candidate_raw)
    last = packet['metadata']['last_mode']
    if (candidate['schema'],candidate['N'],candidate['M'],candidate['last_mode'],
            candidate['target_packet_sha256'],candidate['pivot_result_sha256']) != (
            'ACTUAL_C27_SUPPORT_DIRECTION_V1',N,M,last,target_hash,hashlib.sha256(pivot_raw).hexdigest()):
        raise ValueError('Require the full dyadic direction of this same saved target/support.')
    vq = [fmpq(int(x),1 << candidate['candidate_dyadic_bits'])
          for x in candidate['candidate_dyadic_integers']]
    norm_v = sum((x*x for x in vq),fmpq(0))
    if len(vq) != N+1 or not norm_v > 0:
        raise ValueError('Require a nonzero complete retained direction.')
    scalar,rows,pilot,far = inputs()
    if candidate['basis'] != scalar['basis'] or not N < last < M:
        raise ValueError('Basis/support mismatch.')
    arithmetic = Arithmetic(27,scalar['grid_bits'])
    data = {n:(interval_ball(row['s'],scalar['grid_bits']),interval_ball(row['a'],scalar['grid_bits']))
            for n,row in rows.items()}
    weights = [arb(fmpq(int(x),1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    if len(weights) != M-N or not all(w > 0 for w in weights):
        raise ValueError('Require all positive archived middle weights.')
    selected = list(range(N+1,last+1))
    count = len(selected)
    provenance = dict(N=N,M=M,last_mode=last,basis=scalar['basis'],
                      target_packet_sha256=target_hash,candidate_sha256=hashlib.sha256(candidate_raw).hexdigest(),
                      input_sha256=packet['metadata']['input_sha256'])
    GE,T,F,Bv,weighted_cv,diagonal,qv,cv_energy,finite_hash = finite_action(
        args,selected,vq,scalar,data,weights,arithmetic,provenance)
    kappa = 2/arithmetic.pi
    gram = gram_from_generators(GE,T,diagonal,kappa,selected)
    Kfinite = gram-F
    floor = fmpq(args.floor)
    if not floor > 0:
        raise ValueError('Positive K floor required.')
    generators = arb_mat([[GE[i,j] for j in range(2)]+[T[i,j] for j in range(2)]
                          for i in range(count)])
    floor_verdict,failure,pivots,order = eliminate([Kfinite[i,i]-arb(floor) for i in range(count)],
                    generators,arb_mat(0,0),-kappa,[n*n for n in selected],True,'max-diagonal')
    print('finite K-floor',floor_verdict,'seconds',time.monotonic()-started,flush=True)
    if floor_verdict != 'positive-pivots':
        raise ArithmeticError('Requested finite K floor not certified; no family exclusion.')

    # Reuse exactly the archived extended-support error, on its recorded grid.
    if old_pivots['grid_bits'] != 512 or old_pivots['far_matrix_error'] != packet['metadata']['far_matrix_error']:
        raise ValueError('Wrong archived extended-support far allowance.')
    epsilon = arb(fmpq(int(old_pivots['far_matrix_error'][1]),1 << old_pivots['grid_bits']))
    alpha,scale,Pq,Hq,_ = reconstruct_factors(dict(far,N=last),
                           dict(scalar,rows=[rows[n] for n in range(last+1)]))
    PEq = fmpq_mat([[Pq[i,n] for n in selected] for i in range(Pq.nrows())])
    pvq = fmpq_mat([[sum((Pq[i,n]*vq[n] for n in range(N+1)),fmpq(0))]
                    for i in range(Pq.nrows())])
    PE,pv,core = arb_mat(PEq),arb_mat(pvq),arb(alpha*scale)*arb_mat(Hq)
    Kmid = Kfinite+PE.transpose()*core*PE
    tmid = weighted_cv-Bv+PE.transpose()*core*pv
    inverse_action = Kmid.solve(tmid,algorithm='precond')
    trial_integers = [trial_integer(inverse_action[j,0],args.trial_bits) for j in range(count)]
    xq = [fmpq(n,1 << args.trial_bits) for n in trial_integers]
    x = arb_mat([[arb(t)] for t in xq])
    norm_u = norm_v+sum((t*t for t in xq),fmpq(0))
    pu = arb_mat(pvq-PEq*fmpq_mat([[t] for t in xq]))
    far_center = (pu.transpose()*core*pu)[0,0]
    far_radius = epsilon*arb(norm_u)
    far_energy = far_center+arb(0,far_radius.upper())
    middle_energy = cv_energy-2*(x.transpose()*weighted_cv)[0,0]+(x.transpose()*gram*x)[0,0]
    combined_q = qv-2*(x.transpose()*Bv)[0,0]+(x.transpose()*F*x)[0,0]
    margin = combined_q-middle_energy-far_energy
    gradient = weighted_cv-Bv-Kfinite*x+PE.transpose()*core*pu
    gradient_radius = epsilon*arb(norm_u).sqrt()
    for j in range(count):
        gradient[j,0] += arb(0,gradient_radius.upper())
    gradient_sq = sum((gradient[j,0]**2 for j in range(count)),arb(0))
    correction = gradient_sq/arb(floor)
    optimal_upper = margin+correction
    verdict = ('all-trials-on-this-support-excluded' if optimal_upper < 0 else
               'this-direction-repaired' if margin > 0 else 'inconclusive')
    result = {
        'schema':'ACTUAL_C27_OPTIMAL_SUPPORT_TEST_V1','N':N,'M':M,'last_mode':last,'basis':scalar['basis'],
        'candidate_sha256':hashlib.sha256(candidate_raw).hexdigest(),'finite_action_sha256':finite_hash,
        'target_packet_sha256':target_hash,'pivot_result_sha256':hashlib.sha256(pivot_raw).hexdigest(),
        'working_bits':args.precision,'grid_bits':args.output_bits,'candidate_norm_sq':str(norm_v),
        'combined_norm_sq':str(norm_u),'K_lower_bound':str(floor),
        'K_floor_finite_rows':{'first':N+1,'last':M,'count':M-N},
        'K_floor_positive_pivot_count':len(pivots),'K_floor_pivot_order':[selected[i] for i in order],
        'K_floor_pivot_intervals':[envelope(p,args.output_bits) for p in pivots],
        'trial_dyadic_bits':args.trial_bits,'trial_dyadic_integers':[str(n) for n in trial_integers],
        'combined_Q':envelope(combined_q,args.output_bits),
        'middle_residual_energy':envelope(middle_energy,args.output_bits),
        'far_residual_energy':envelope(far_energy,args.output_bits),
        'archived_far_error_upper':envelope(epsilon,args.output_bits),
        'gradient_intervals':[envelope(gradient[j,0],args.output_bits) for j in range(count)],
        'fixed_trial_margin':envelope(margin,args.output_bits),
        'optimizer_correction_upper':envelope(correction,args.output_bits),
        'optimized_family_quadratic_upper':envelope(optimal_upper,args.output_bits),
        'verdict':verdict,'elapsed_seconds':time.monotonic()-started,
        'input_sha256':packet['metadata']['input_sha256'],
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'python_flint':flint.__version__,'FLINT':flint.__FLINT_VERSION__,
        'conditions':['Same actual complete Fourier/common-form and exterior-domain identifications.',
                      'Same positive actual exterior H and archived W majorizing H^-1.',
                      'Sound reused target/pivot data and extended-support far error.'],
        'scope':'The existing square-completion bound m(x)+||Tv-Kx||²/kappa bounds the best lower target '
                'on this same full retained direction among ALL trials supported in modes1281..last_mode with this W. '
                'A strictly negative upper excludes that whole sufficient family, not actual Weil positivity or other supports/majorants. '
                'A positive fixed margin repairs only this direction. No old scientific producer reruns, generic algorithm, Lean, cofinal or RH result.',
    }
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(verdict,'margin',margin.str(24),'correction',correction.str(16),
          'optimized upper',optimal_upper.str(24),'seconds',time.monotonic()-started,flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare',action='store_true')
    parser.add_argument('--self-check',action='store_true')
    parser.add_argument('--target',type=Path,default=HERE/'structured-joint-2048-target.json.zlib')
    parser.add_argument('--pivots',type=Path,default=HERE/'structured-joint-2048-pivots.json')
    parser.add_argument('--candidate',type=Path)
    parser.add_argument('--finite',type=Path,help='Reuse an existing finite Gram/action packet; do not rescan middle rows.')
    parser.add_argument('--finite-out',type=Path,help='Archive a new finite Gram/action packet before factorization.')
    parser.add_argument('--mode',type=int,default=0)
    parser.add_argument('--floor',default='1/16')
    parser.add_argument('--precision',type=int,default=768)
    parser.add_argument('--output-bits',type=int,default=512)
    parser.add_argument('--direction-bits',type=int,default=256)
    parser.add_argument('--trial-bits',type=int,default=384)
    parser.add_argument('--out',type=Path)
    args = parser.parse_args()
    if args.precision < max(640,args.output_bits,args.direction_bits,args.trial_bits)+96:
        parser.error('Use at least96 guard bits over source/output/dyadic grids.')
    if (args.finite and args.finite_out) or (args.finite_out and args.finite_out.exists()):
        parser.error('Finite input/output must be distinct and output fresh.')
    ctx.prec = args.precision
    if args.self_check:
        self_check()
    elif args.out is None or args.out.exists() or (not args.prepare and args.candidate is None):
        parser.error('Use a fresh output and a saved direction unless preparing it.')
    elif args.prepare:
        prepare(args)
    else:
        certify(args)


if __name__ == '__main__':
    main()
