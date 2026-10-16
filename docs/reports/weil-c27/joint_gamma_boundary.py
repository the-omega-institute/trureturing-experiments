"""Apply the actual positive Gamma channels to every retained coefficient.

Reuse the saved complete residual target and its exact two-generator trial.
Only previously unsaved exterior generators and thin boundary projections
are produced. The correlated far Gram is reused with a matrix Young bound.
The existing signed-displacement elimination consumes the new lower target.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import zlib

import flint
from flint import arb, arb_mat, arb_poly, ctx, fmpq

from arb_arithmetic import Arithmetic
from arb_inverse_residual_obstruction import HERE, M, N, inputs
from diagonal_floor_witness import interval_ball
from gamma_boundary_update import boundary_columns, fibonacci_bins, read_json
from structured_inverse_residual_interval import eliminate, envelope
from structured_residual_schur_finish import load_target


def matrix_in(rows, bits):
    return arb_mat([[interval_ball(x,bits) for x in row] for row in rows])


def matrix_out(A, bits):
    return [[envelope(A[i,j],bits) for j in range(A.ncols())] for i in range(A.nrows())]


def write_packet(path, packet):
    raw=zlib.compress(json.dumps(packet,separators=(',',':')).encode(),9)
    path.write_bytes(raw)
    return hashlib.sha256(raw).hexdigest()


def original_target(args):
    packet,_,_,target_hash,pivot_raw=load_target(args.target,args.pivots)
    if packet['metadata']['last_mode'] != 2048:
        raise ValueError('This consumer uses the existing exact trial through2048.')
    proof=json.loads(pivot_raw)
    for key in ('far_matrix_error','complete_far_operator_allowance','trial_operator_norm_sq_upper'):
        if packet['metadata'][key] != proof[key]:
            raise ValueError('Original metadata and its explicitly gridded result disagree.')
    return packet,target_hash,proof['grid_bits']


def generator_action(args, packet, target_hash, scalar, rows, arithmetic):
    """New unsaved GO=G-H_E F^-1 G_E, with existing geometric error bounds."""
    last=packet['metadata']['last_mode']
    provenance=dict(target_packet_sha256=target_hash,input_sha256=packet['metadata']['input_sha256'],
                    N=N,M=M,last_mode=last,basis=packet['metadata']['basis'])
    if args.generator:
        saved,gh=read_json(args.generator)
        if saved['schema'] != 'ACTUAL_C27_RESIDUAL_GENERATOR_V1' or saved['metadata'] != provenance:
            raise ValueError('Wrong actual residual-generator packet.')
        GO=matrix_in(saved['GO'],saved['grid_bits'])
        if (GO.nrows(),GO.ncols()) != (M-last,2):
            raise ValueError('Incomplete actual exterior generators.')
        return GO,gh
    started=time.monotonic()
    finite,fh=read_json(args.support)
    if finite['metadata']['target_packet_sha256'] != target_hash or finite['metadata']['last_mode'] != last:
        raise ValueError('Finite selected block belongs to another actual target.')
    bits=finite['grid_bits']
    GE=matrix_in(finite['GE'],bits)
    selected=list(range(N+1,last+1))
    kap=interval_ball(packet['kappa'],packet['grid_bits'])
    F=arb_mat(len(selected),len(selected))
    for i,m in enumerate(selected):
        F[i,i]=interval_ball(finite['F_diagonal'][i],bits)
        for j,n in enumerate(selected[:i]):
            F[i,j]=kap*(GE[i,0]*GE[j,1]-GE[i,1]*GE[j,0])/(m*m-n*n)
            F[j,i]=F[i,j]
    inverse=F.solve(GE,algorithm='precond')
    data={n:(interval_ball(r['s'],scalar['grid_bits']),interval_ball(r['a'],scalar['grid_bits']))
          for n,r in rows.items()}
    a=[[GE[i,1]*inverse[i,k] for i in range(len(selected))] for k in range(2)]
    b=[[GE[i,0]*inverse[i,k] for i in range(len(selected))] for k in range(2)]
    sums_a=[sum(map(abs,v),arb(0)) for v in a]
    sums_b=[sum(map(abs,v),arb(0)) for v in b]
    order=384
    ratios=[arb(fmpq(n*n,last*last)) for n in selected]
    power=[arb(1)]*len(selected)
    A,B=[[] for _ in range(2)],[[] for _ in range(2)]
    for _ in range(order):
        for k in range(2):
            A[k].append(sum((a[k][i]*power[i] for i in range(len(selected))),arb(0)))
            B[k].append(sum((b[k][i]*power[i] for i in range(len(selected))),arb(0)))
        power=[power[i]*ratios[i] for i in range(len(selected))]
    pa,pb=list(map(arb_poly,A)),list(map(arb_poly,B))
    GO=arb_mat(M-last,2)
    checks=[]
    for m in range(last+1,M+1):
        gm0,gm1=arb((-1)**m),arb((-1)**m)*m*data[m][0]
        if m <= 2*last:
            z=[arithmetic.even_entry(m,n,data) for n in selected]
            action=[sum((z[i]*inverse[i,k] for i in range(len(selected))),arb(0)) for k in range(2)]
        else:
            t=arb(fmpq(last*last,m*m))
            action=[]
            for k in range(2):
                value=kap/(m*m)*(gm0*pa[k](t)-gm1*pb[k](t))
                error=kap/(m*m)*t**order/(1-t)*(abs(gm0)*sums_a[k]+abs(gm1)*sums_b[k])
                action.append(value+arb(0,error.upper()))
            if m in (2*last+1,8192,16384,M):
                z=[arithmetic.even_entry(m,n,data) for n in selected]
                direct=[sum((z[i]*inverse[i,k] for i in range(len(selected))),arb(0)) for k in range(2)]
                if not all(action[k].overlaps(direct[k]) for k in range(2)):
                    raise ArithmeticError('New generator action disagrees with direct actual rows.')
                checks.append(m)
        GO[m-last-1,0],GO[m-last-1,1]=gm0-action[0],gm1-action[1]
    saved=dict(schema='ACTUAL_C27_RESIDUAL_GENERATOR_V1',metadata=provenance,
               grid_bits=args.output_bits,working_bits=args.precision,
               GO=matrix_out(GO,args.output_bits),selected_inverse_generator=matrix_out(inverse,args.output_bits),
               support_action_sha256=fh,direct_last_mode=2*last,geometric_order=order,
               direct_sample_checks=checks,elapsed_seconds=time.monotonic()-started,
               source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if args.generator_out is None:
        raise ValueError('Archive the new residual generators for future consumers.')
    gh=write_packet(args.generator_out,saved)
    print('new residual generators archived; seconds',time.monotonic()-started,flush=True)
    return matrix_in(saved['GO'],args.output_bits),gh


def produce(args):
    started=time.monotonic()
    packet,th,_=original_target(args)
    scalar,rows,pilot,_=inputs()
    arithmetic=Arithmetic(27,scalar['grid_bits'])
    last=packet['metadata']['last_mode']
    GO,gh=generator_action(args,packet,th,scalar,rows,arithmetic)
    weights=[arb(fmpq(int(x),1<<pilot['middle_weight_bits'])) for x in pilot['middle_weight_numerators']]
    D=arb(fmpq(int(pilot['far_weight_numerator']),int(pilot['far_weight_denominator'])))
    bins=fibonacci_bins(args.last_bin)
    U,eta=boundary_columns(arithmetic,bins,D)
    WU=arb_mat([[weights[i]*U[i,k] for k in range(len(bins))] for i in range(M-N)])
    A=U.transpose()*WU
    order=512
    moments=arb_mat(order,2*len(bins))
    absolute=arb_mat(1,2*len(bins))
    for first in range(last+1,M+1,256):
        modes=list(range(first,min(first+256,M+1)))
        powers=arb_mat(len(modes),order)
        V=arb_mat(len(modes),2*len(bins))
        for i,m in enumerate(modes):
            ratio=arb(fmpq(N*N,m*m))
            value=arb(1)
            for r in range(order):
                powers[i,r]=value
                value *= ratio
            for k in range(len(bins)):
                for j in range(2):
                    v=WU[m-N-1,k]*GO[m-last-1,j]/(m*m)
                    V[i,2*k+j]=v
                    absolute[0,2*k+j] += abs(v)
        moments += powers.transpose()*V
        if (first-last-1)%4096 == 0:
            print('new thin finite moments through',modes[-1],'seconds',time.monotonic()-started,flush=True)
    GP=arb_mat([[interval_ball(x,packet['grid_bits']) for x in row[:2]] for row in packet['generators']])
    kap=interval_ball(packet['kappa'],packet['grid_bits'])
    q=arb(fmpq(N,last+1))
    geometric=q**(2*order)/(1-q*q)
    beta=arb_mat(len(bins),N+1)
    for k in range(len(bins)):
        a=arb_poly([moments[r,2*k] for r in range(order)])
        b=arb_poly([moments[r,2*k+1] for r in range(order)])
        for n in range(N+1):
            t=arb(fmpq(n*n,N*N))
            v=kap*(GP[n,1]*a(t)-GP[n,0]*b(t))
            error=kap*geometric*(abs(GP[n,1])*absolute[0,2*k]+abs(GP[n,0])*absolute[0,2*k+1])
            beta[k,n]=v+arb(0,error.upper())
    result=dict(schema='ACTUAL_C27_JOINT_GAMMA_ACTION_V1',N=N,M=M,last_mode=last,
                target_packet_sha256=th,generator_packet_sha256=gh,basis=packet['metadata']['basis'],
                grid_bits=args.output_bits,working_bits=args.precision,bins=bins,
                finite_gram=matrix_out(A,args.output_bits),finite_beta=matrix_out(beta,args.output_bits),
                far_gram_trace_bound=envelope(sum(eta,arb(0)),args.output_bits),
                finite_moment_order=order,finite_moment_ratio=envelope(q,args.output_bits),
                finite_geometric_factor=envelope(geometric,args.output_bits),
                finite_absolute_moments=matrix_out(absolute,args.output_bits),
                input_sha256=packet['metadata']['input_sha256'],
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                boundary_source_sha256=hashlib.sha256((HERE/'gamma_boundary_update.py').read_bytes()).hexdigest(),
                elapsed_seconds=time.monotonic()-started)
    write_packet(args.out,result)
    print('new full-coefficient thin action archived; seconds',time.monotonic()-started,flush=True)


def consume(args):
    started=time.monotonic()
    packet,th,metadata_bits=original_target(args)
    finite,fh=read_json(args.finite)
    if (finite['schema'],finite['target_packet_sha256'],finite['N'],finite['M'],finite['last_mode']) != (
            'ACTUAL_C27_JOINT_GAMMA_ACTION_V1',th,N,M,packet['metadata']['last_mode']):
        raise ValueError('Wrong complete finite boundary projections.')
    rho=fmpq(args.rho)
    if not rho > 0:
        raise ValueError('The common Young parameter must be positive.')
    bits=finite['grid_bits']
    beta,A=matrix_in(finite['finite_beta'],bits),matrix_in(finite['finite_gram'],bits)
    count=len(finite['bins'])
    if (beta.nrows(),beta.ncols(),A.nrows(),A.ncols()) != (count,N+1,count,count):
        raise ValueError('Incomplete common coefficient map.')
    delta=interval_ball(finite['far_gram_trace_bound'],bits).upper()
    S=arb_mat(A)
    for k in range(count):
        S[k,k] += 1+(1+1/arb(rho))*delta
    inverse=S.inv()
    action=inverse*beta
    gain=beta.transpose()*action
    oldbits=packet['grid_bits']
    oldG=matrix_in(packet['generators'],oldbits)
    oldK=matrix_in(packet['kernel'],oldbits)
    oldrank=oldK.nrows()
    P=arb_mat([[oldG[n,4+k] for n in range(N+1)] for k in range(oldrank)])
    far=P.transpose()*oldK*P
    allowance=interval_ball(packet['metadata']['complete_far_operator_allowance'],metadata_bits).upper()
    diagonal=[interval_ball(x,oldbits)+gain[n,n]-arb(rho)*(far[n,n]+allowance)
              for n,x in enumerate(packet['diagonal'])]
    kernel=arb_mat(oldrank+count,oldrank+count)
    for i in range(oldrank):
        for j in range(oldrank):
            kernel[i,j]=(1+arb(rho))*oldK[i,j]
    for i in range(count):
        for j in range(count):
            kernel[oldrank+i,oldrank+j]=-inverse[i,j]
    features=arb_mat([[oldG[n,4+k] for k in range(oldrank)]+[beta[k,n] for k in range(count)]
                      for n in range(N+1)])
    generators=arb_mat([[oldG[n,k] for k in range(4)]+[features[n,k] for k in range(oldrank+count)]
                       +[n*n*features[n,k] for k in range(oldrank+count)] for n in range(N+1)])
    metadata=dict(packet['metadata'],joint_boundary_action_sha256=fh,original_target_sha256=th,
                  positive_boundary_bins=finite['bins'],common_young_parameter=str(rho),
                  original_metadata_grid_bits=metadata_bits,
                  reused_far_allowance=envelope(allowance,args.output_bits),
                  reused_far_allowance_grid_bits=args.output_bits,
                  signed_update_rank=oldrank+count,total_generator_columns=generators.ncols())
    updated=dict(schema='ACTUAL_C27_STRUCTURED_TARGET_V1',grid_bits=args.output_bits,
                 metadata=metadata,producer_working_bits=args.precision,
                 producer_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 diagonal=[envelope(x,args.output_bits) for x in diagonal],
                 generators=matrix_out(generators,args.output_bits),kernel=matrix_out(kernel,args.output_bits),
                 kappa=packet['kappa'])
    # The original kappa pair uses the original grid; transport it if needed.
    kap=interval_ball(packet['kappa'],oldbits)
    updated['kappa']=envelope(kap,args.output_bits)
    if args.target_out is None:
        raise ValueError('Archive the new target before elimination; do not repeat thin production.')
    uh=write_packet(args.target_out,updated)
    # Apply the certificate to the same outward target saved for later consumers.
    diagonal=[interval_ball(x,args.output_bits) for x in updated['diagonal']]
    generators=matrix_in(updated['generators'],args.output_bits)
    kernel=matrix_in(updated['kernel'],args.output_bits)
    kap=interval_ball(updated['kappa'],args.output_bits)
    verdict,failure,pivots,order=eliminate(diagonal,generators,kernel,kap,[n*n for n in range(N+1)],True,'max-diagonal')
    result=dict(schema='ACTUAL_C27_STRUCTURED_JOINT_INTERVAL_V1',N=N,M=M,last_mode=packet['metadata']['last_mode'],
                basis=packet['metadata']['basis'],target_packet_sha256=uh,original_target_packet_sha256=th,
                finite_boundary_action_sha256=fh,common_young_parameter=str(rho),
                working_bits=args.precision,grid_bits=args.output_bits,verdict=verdict,
                failure_mode=failure,pivot_count=len(pivots),pivot_modes=order,
                pivot_intervals=[envelope(p,args.output_bits) for p in pivots],
                elapsed_seconds=time.monotonic()-started,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                scope='The same exact finite inverse trial and actual exterior, with nineteen positive Gamma boundary channels '
                      'and a common correlated far Young bound. Positive pivots would certify this full sufficient lower target; '
                      'negative pivots reject this target only. No old scalar, Gram, target or moment producer rerun; '
                      'no generic theorem, Lean, cofinal positivity or RH claim.')
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(verdict,'positive count',len(pivots)-(verdict!='positive-pivots'),'failure',failure,
          'seconds',time.monotonic()-started,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--produce',action='store_true')
    p.add_argument('--target',type=Path,default=HERE/'structured-joint-2048-target.json.zlib')
    p.add_argument('--pivots',type=Path,default=HERE/'structured-joint-2048-pivots.json')
    p.add_argument('--support',type=Path,default=HERE/'support-2048-finite-action.json.zlib')
    p.add_argument('--generator',type=Path)
    p.add_argument('--generator-out',type=Path)
    p.add_argument('--finite',type=Path)
    p.add_argument('--last-bin',type=int,default=6765)
    p.add_argument('--rho',default='1/32')
    p.add_argument('--precision',type=int,default=1024)
    p.add_argument('--output-bits',type=int,default=640)
    p.add_argument('--target-out',type=Path)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    outputs=[x for x in [a.out,a.generator_out,a.target_out] if x is not None]
    if any(x.exists() for x in outputs) or len(set(outputs)) != len(outputs):
        p.error('Every output must be fresh and distinct.')
    if a.precision < max(640,a.output_bits)+192 or (not a.produce and a.finite is None):
        p.error('Use at least192 guard bits and a saved finite action unless producing it.')
    if a.generator is not None and a.generator_out is not None:
        p.error('Reuse existing generators or produce new ones, not both.')
    ctx.prec=a.precision
    produce(a) if a.produce else consume(a)


if __name__ == '__main__':
    main()
