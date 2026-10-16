"""Enclose a full c27 inverse-residual target through finite generators.

The only solve has two RHS. The selected inverse trial is exact F^-1 B,
not the rounded floating proposal. All finite middle rows and the existing
infinite far moments enter the same coefficient-uniform lower target.
Finite displacement elimination reports positive, negative or inconclusive
pivots; none is a Lean or cofinal RH result.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
from types import SimpleNamespace
import zlib

import flint
from flint import arb, arb_mat, ctx, fmpq

from arb_arithmetic import Arithmetic
from arb_far_coupling import far_bounds, reconstruct_factors
from arb_inverse_residual_obstruction import HERE, M, N, inputs
from diagonal_floor_witness import interval_ball


def envelope(value, bits):
    """Report actual outward width, including an inconclusive wide interval."""
    if not value.is_finite():
        return None
    scale = 1 << bits
    return [str(int((value.lower()*scale).floor().unique_fmpz())),
            str(int((value.upper()*scale).ceil().unique_fmpz()))]


def eliminate(diagonal, generators, kernel, kappa, nodes, progress=False, pivoting='natural'):
    """Finite symmetric displacement LDL; stop before any uncertain division.

    Generator columns are [GP,T,P,D*P]. The skew form is
    kappa*[[J,-J],[-J,0]] plus [[0,kernel],[-kernel,0]].
    Diagonals remain independently supplied and updated.
    """
    size, width = generators.nrows(), generators.ncols()
    rank = (width-4)//2
    if width != 4+2*rank or len(diagonal) != size or len(nodes) != size:
        raise ValueError('Invalid finite generator shape.')
    pivots = []
    order = []
    work = list(diagonal)
    G = generators
    remaining = list(range(size))
    for step in range(size):
        k = remaining[0] if pivoting == 'natural' else max(remaining,key=lambda i:float(work[i].mid()))
        pivot = work[k]
        pivots.append(pivot)
        order.append(k)
        if not pivot > 0:
            return ('negative-pivot' if pivot < 0 else 'inconclusive-pivot'), k, pivots, order
        remaining.remove(k)
        if not remaining:
            break
        row = [G[k,j] for j in range(width)]
        action = [kappa*(row[1]-row[3]), -kappa*(row[0]-row[2]),
                  -kappa*row[1], kappa*row[0]]
        if rank:
            p = arb_mat([[row[4+j]] for j in range(rank)])
            dp = arb_mat([[row[4+rank+j]] for j in range(rank)])
            kp, kdp = kernel*p, kernel*dp
            action += [kdp[j,0] for j in range(rank)]
            action += [-kp[j,0] for j in range(rank)]
        raw = G*arb_mat([[x] for x in action])
        multipliers = [arb(0)]*size
        for i in remaining:
            entry = raw[i,0]/(nodes[i]-nodes[k])
            multipliers[i] = entry/pivot
            work[i] -= entry*multipliers[i]
        G -= arb_mat([[x] for x in multipliers])*arb_mat([row])
        if progress and (step % 64 == 0):
            print('positive pivot count',step+1,'mode',k,'pivot',pivot.str(12),flush=True)
    return 'positive-pivots', None, pivots, order


def self_check():
    from flint import fmpq_mat
    J = fmpq_mat([[0,1],[-1,0]])
    selected, retained, outside = [5,7,9], [0,1,2,3], [11,12,13,14]
    kap = fmpq(2,3)
    def g(n):
        return [fmpq((-1)**n), fmpq((-1)**n*(n+1),n+2)]
    def block(ms,ns):
        return fmpq_mat([[fmpq(50+n) if m == n else
            kap*(g(m)[0]*g(n)[1]-g(m)[1]*g(n)[0])/(m*m-n*n)
            for n in ns] for m in ms])
    def gen(ms):
        return fmpq_mat([g(m) for m in ms])
    def diag(xs):
        return fmpq_mat([[x if i == j else 0 for j in range(len(xs))] for i,x in enumerate(xs)])
    F, B, Z = block(selected,selected), block(selected,retained), block(outside,selected)
    U, Y = F.solve(gen(selected)), F.solve(B)
    GP, GO = gen(retained)-B.transpose()*U, gen(outside)-Z*U
    reconstructed = fmpq_mat([[kap*(U[i,0]*GP[j,1]-U[i,1]*GP[j,0])/(m*m-n*n)
                              for j,n in enumerate(retained)] for i,m in enumerate(selected)])
    assert reconstructed == Y
    schur = block(retained,retained)-B.transpose()*Y
    R = block(outside,retained)-Z*Y
    reconstructed_R = fmpq_mat([[kap*(GO[i,0]*GP[j,1]-GO[i,1]*GP[j,0])/(m*m-n*n)
                                for j,n in enumerate(retained)] for i,m in enumerate(outside)])
    assert reconstructed_R == R
    W = diag([fmpq(1,i+1) for i in range(len(outside))])
    gram, T = R.transpose()*W*R, R.transpose()*W*GO
    D = diag([n*n for n in retained])
    assert D*schur-schur*D == kap*GP*J*GP.transpose()
    assert D*gram-gram*D == kap*(T*J*GP.transpose()+GP*J*T.transpose())
    target = schur-gram
    P = fmpq_mat([[fmpq(n+1,20),fmpq(n+2,30)] for n in retained])
    K = fmpq_mat([[1,fmpq(1,10)],[fmpq(1,10),fmpq(-1,2)]])
    full = target-P*K*P.transpose()
    G = arb_mat([[arb(GP[i,j]) for j in range(2)]+[arb(T[i,j]) for j in range(2)]
                 +[arb(P[i,j]) for j in range(2)]+[arb((D*P)[i,j]) for j in range(2)]
                 for i in range(len(retained))])
    verdict, failure, pivots, order = eliminate([arb(full[i,i]) for i in range(len(retained))],
                                        G,arb_mat(K),arb(kap),[n*n for n in retained])
    assert verdict == 'positive-pivots' and failure is None
    # Independently computed exact dense LDL pivots, including a signed update.
    dense = fmpq_mat(full)
    for k in range(len(retained)):
        exact = dense[k,k]
        assert pivots[k].contains(arb(exact))
        for i in range(k+1,len(retained)):
            for j in range(k+1,len(retained)):
                dense[i,j] -= dense[i,k]*dense[k,j]/exact
    pivoted = eliminate([arb(full[i,i]) for i in range(len(retained))],
                         G,arb_mat(K),arb(kap),[n*n for n in retained],pivoting='max-diagonal')
    assert pivoted[0] == 'positive-pivots' and sorted(pivoted[3]) == list(range(len(retained)))
    negative = [arb(full[i,i]) for i in range(len(retained))]
    negative[0] = arb(-1)
    assert eliminate(negative,G,arb_mat(K),arb(kap),[n*n for n in retained])[0] == 'negative-pivot'
    print('exact rational Schur/residual/Gram and dense-LDL comparison passed; negative diagonal rejected',flush=True)


def build_target(args):
    started = time.monotonic()
    scalar, rows, pilot, far = inputs()
    last = args.last_mode
    if not N < last < M:
        raise ValueError('Selected modes must be in the published finite middle.')
    arithmetic = Arithmetic(27, scalar['grid_bits'])
    kappa = 2/arithmetic.pi
    data = {n:(interval_ball(r['s'],scalar['grid_bits']),interval_ball(r['a'],scalar['grid_bits']))
            for n,r in rows.items()}
    modes = list(range(last+1))
    selected = list(range(N+1,last+1))
    def g(n):
        phase = arb(-1 if n % 2 else 1)
        if n == 0:
            phase /= arb(2).sqrt()
        return [phase,phase*n*data[n][0]]
    original = {n:g(n) for n in range(M+1)}
    F = arb_mat([[arithmetic.even_entry(m,n,data) for n in selected] for m in selected])
    B = arb_mat([[arithmetic.even_entry(m,n,data) for n in range(N+1)] for m in selected])
    U = F.solve(arb_mat([original[n] for n in selected]))
    transformed = arb_mat([original[n] for n in range(N+1)])-B.transpose()*U
    GP = [[transformed[i,j] for j in range(2)] for i in range(N+1)]
    Y = arb_mat([[kappa*(U[i,0]*GP[n][1]-U[i,1]*GP[n][0])/(m*m-n*n)
                  for n in range(N+1)] for i,m in enumerate(selected)])
    schur_diag = [arithmetic.even_entry(n,n,data)
                  -sum((B[j,n]*Y[j,n] for j in range(len(selected))),arb(0)) for n in range(N+1)]
    print('directed two-RHS trial reconstructed; seconds',time.monotonic()-started,flush=True)
    weights = [arb(fmpq(int(x),1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    T = arb_mat(N+1,2)
    gram_diag = [arb(0)]*(N+1)
    for first in range(last+1,M+1,128):
        outside = list(range(first,min(first+128,M+1)))
        Z = arb_mat([[arithmetic.even_entry(m,n,data) for n in selected] for m in outside])
        GO = arb_mat([original[m] for m in outside])-Z*U
        residual = arb_mat([[kappa*(GO[i,0]*GP[n][1]-GO[i,1]*GP[n][0])/(m*m-n*n)
                             for n in range(N+1)] for i,m in enumerate(outside)])
        WG = arb_mat([[weights[m-N-1]*GO[i,j] for j in range(2)] for i,m in enumerate(outside)])
        T += residual.transpose()*WG
        for n in range(N+1):
            gram_diag[n] += sum((weights[m-N-1]*residual[i,n]**2 for i,m in enumerate(outside)),arb(0))
        if (first-last-1) % 4096 == 0:
            print('directed finite middle through',outside[-1],'seconds',time.monotonic()-started,flush=True)

    alpha, scale, Pq, Hq, _ = reconstruct_factors(
        dict(far,N=last),dict(scalar,rows=[rows[n] for n in modes]))
    feature = arb_mat([[arb(Pq[i,n]) for n in range(N+1)] for i in range(Pq.nrows())])
    feature -= arb_mat([[arb(Pq[i,n]) for n in selected] for i in range(Pq.nrows())])*Y
    K = arb(alpha*scale)*arb_mat(Hq)
    KF = K*feature
    far_diag = [sum((feature[j,n]*KF[j,n] for j in range(feature.nrows())),arb(0))
                for n in range(N+1)]
    symbol_bound = 2*arithmetic.H+sum((w for w,_ in arithmetic.primes),arb(0))
    scalar_radius = max(arb(int(rows[n]['s'][1])-int(rows[n]['s'][0]))
                        /(1 << (scalar['grid_bits']+1)) for n in modes)
    bounds = far_bounds(SimpleNamespace(M=M,L=arithmetic.d,pi=arithmetic.pi,budget=symbol_bound),
                        last,far['jet_order'],scalar_radius,far['grid_bits'],far['moment_diagonals'])
    epsilon = arb(scale)*bounds['unweighted_matrix_error']
    trial_norm_bound = 1+sum((v**2 for v in Y.entries()),arb(0))
    allowance = arb((epsilon*trial_norm_bound).upper())
    full_diagonal = [schur_diag[n]-gram_diag[n]-far_diag[n]-allowance for n in range(N+1)]
    rank = feature.nrows()
    generators = arb_mat([GP[n]+[T[n,j] for j in range(2)]
                          +[feature[j,n] for j in range(rank)]
                          +[n*n*feature[j,n] for j in range(rank)] for n in range(N+1)])
    print('complete target generators enclosed; seconds',time.monotonic()-started,flush=True)
    metadata = {
        'N':N,'M':M,'last_mode':last,'basis':scalar['basis'],
        'finite_inverse_rhs_count':2,'selected_exact_residual_zero_rows':len(selected),
        'finite_nonzero_residual_rows':{'first':last+1,'last':M,'count':M-last},
        'far_feature_rank':rank,'total_generator_columns':generators.ncols(),
        'far_matrix_error':envelope(epsilon,args.output_bits),
        'trial_operator_norm_sq_upper':envelope(trial_norm_bound,args.output_bits),
        'complete_far_operator_allowance':envelope(allowance,args.output_bits),
        'input_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['arithmetic-scalars.json.zlib','middle-scalars.json.zlib','diagonal-floor-witness.json',
             'far-coupling.json.zlib','inverse-residual-pilot.json.zlib']},
    }
    return full_diagonal,generators,K,kappa,metadata


def run(args):
    started = time.monotonic()
    packet_hash = None
    if args.target:
        raw = args.target.read_bytes()
        packet_hash = hashlib.sha256(raw).hexdigest()
        packet = json.loads(zlib.decompress(raw))
        if (packet['schema'],packet['metadata']['N'],packet['metadata']['M'],packet['metadata']['last_mode']) != (
                'ACTUAL_C27_STRUCTURED_TARGET_V1',N,M,args.last_mode):
            raise ValueError('Wrong full-target packet.')
        for name,digest in packet['metadata']['input_sha256'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Target input pin mismatch: '+name)
        bits = packet['grid_bits']
        def decode_matrix(rows):
            return arb_mat([[interval_ball(x,bits) for x in row] for row in rows])
        full_diagonal = [interval_ball(x,bits) for x in packet['diagonal']]
        generators = decode_matrix(packet['generators'])
        K = decode_matrix(packet['kernel'])
        kappa = interval_ball(packet['kappa'],bits)
        metadata = packet['metadata']
        if (len(full_diagonal),generators.nrows(),generators.ncols(),K.nrows(),K.ncols()) != (
                N+1,N+1,metadata['total_generator_columns'],metadata['far_feature_rank'],metadata['far_feature_rank']):
            raise ValueError('Incomplete target dimensions.')
        print('reused complete directed target; no middle/far production rerun',flush=True)
    else:
        full_diagonal,generators,K,kappa,metadata = build_target(args)
        if args.target_out:
            packet_bits = 640
            packet = {
                'schema':'ACTUAL_C27_STRUCTURED_TARGET_V1','grid_bits':packet_bits,
                'metadata':metadata,'producer_working_bits':args.precision,
                'producer_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'diagonal':[envelope(x,packet_bits) for x in full_diagonal],
                'generators':[[envelope(generators[i,j],packet_bits) for j in range(generators.ncols())]
                              for i in range(generators.nrows())],
                'kernel':[[envelope(K[i,j],packet_bits) for j in range(K.ncols())] for i in range(K.nrows())],
                'kappa':envelope(kappa,packet_bits),
            }
            raw = zlib.compress(json.dumps(packet,separators=(',',':')).encode(),9)
            args.target_out.write_bytes(raw)
            packet_hash = hashlib.sha256(raw).hexdigest()
            # Judge the same outward exported target that later consumers reuse.
            full_diagonal = [interval_ball(x,packet_bits) for x in packet['diagonal']]
            generators = arb_mat([[interval_ball(x,packet_bits) for x in row] for row in packet['generators']])
            K = arb_mat([[interval_ball(x,packet_bits) for x in row] for row in packet['kernel']])
            kappa = interval_ball(packet['kappa'],packet_bits)
            print('reusable complete target exported; bytes',len(raw),'seconds',time.monotonic()-started,flush=True)
    verdict,failure,pivots,order = eliminate(full_diagonal,generators,K,kappa,
                                            [n*n for n in range(N+1)],True,args.pivoting)
    result = dict(metadata,**{
        'schema':'ACTUAL_C27_STRUCTURED_JOINT_INTERVAL_V1',
        'working_bits':args.precision,'grid_bits':args.output_bits,
        'verdict':verdict,'first_nonpositive_or_uncertain_pivot':failure,
        'pivot_intervals':[envelope(p,args.output_bits) for p in pivots],
        'pivot_modes':order,'pivoting':args.pivoting,'target_packet_sha256':packet_hash,
        'elapsed_seconds':time.monotonic()-started,
        'python_flint':flint.__version__,'FLINT':flint.__FLINT_VERSION__,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'conditions':['Existing complete actual Fourier/common-form and finite exterior operator-domain premises.',
                      'Positive actual exterior H with the same archived W majorizing H^-1.'],
        'scope':'Coefficient-uniform sufficient lower target for the actual c27 Schur form, using exact selected finite inverse F^-1 B. '
                'Positive pivots would certify this numerical lower target under the stated premises; a negative pivot would reject this target, not the actual Weil form. '
                'An inconclusive pivot certifies neither sign. Actual diagonals, all finite middle residual rows, the existing far core and its paid operator allowance are included. '
                'Finite generators only; no bounded infinite displacement claim, new general theorem, Lean verification, cofinal positivity or RH result.',
    })
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(verdict,'pivot',failure,'seconds',time.monotonic()-started,flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--last-mode',type=int,default=2048)
    parser.add_argument('--precision',type=int,default=768)
    parser.add_argument('--output-bits',type=int,default=512)
    parser.add_argument('--self-check',action='store_true')
    parser.add_argument('--pivoting',choices=['natural','max-diagonal'],default='max-diagonal')
    parser.add_argument('--target',type=Path,help='Reuse this complete outward target; run only the factorization.')
    parser.add_argument('--target-out',type=Path,help='Export complete target for subsequent consumers.')
    parser.add_argument('--out',type=Path)
    args = parser.parse_args()
    if args.precision < max(640,args.output_bits)+96:
        parser.error('Keep at least96 guard bits over the source/output grids.')
    ctx.prec = args.precision
    if args.self_check:
        self_check()
        return
    if args.out is None or args.out.exists():
        parser.error('Use a fresh output path.')
    if (args.target and args.target_out) or (args.target_out and args.target_out.exists()):
        parser.error('Target input/output must be distinct and output fresh.')
    run(args)


if __name__ == '__main__':
    main()
