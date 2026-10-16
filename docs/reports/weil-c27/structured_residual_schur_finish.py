"""Reuse a complete directed target and solve its unresolved Schur block.

The positive coordinate block is identified by the saved interval pivot
certificate. FLINT's validated dense solve and small-matrix eigensolver
replace further generator elimination. No target producer is rerun.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import zlib

import flint
from flint import arb, arb_mat, ctx

from arb_inverse_residual_obstruction import HERE, M, N
from diagonal_floor_witness import interval_ball
from structured_inverse_residual_interval import envelope


def run(args):
    started = time.monotonic()
    raw = args.target.read_bytes()
    target_hash = hashlib.sha256(raw).hexdigest()
    packet = json.loads(zlib.decompress(raw))
    pivot_raw = args.pivots.read_bytes()
    proof = json.loads(pivot_raw)
    if (packet['schema'],packet['metadata']['N'],packet['metadata']['M'],proof['schema'],
            proof['target_packet_sha256'],proof['verdict']) != (
            'ACTUAL_C27_STRUCTURED_TARGET_V1',N,M,'ACTUAL_C27_STRUCTURED_JOINT_INTERVAL_V1',
            target_hash,'inconclusive-pivot'):
        raise ValueError('Require the same complete target and its incomplete pivot certificate.')
    for name,digest in packet['metadata']['input_sha256'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Actual input pin mismatch: '+name)
    prefix = proof['pivot_modes'][:-1]
    if len(set(prefix)) != len(prefix) or not all(0 <= n <= N for n in prefix):
        raise ValueError('Invalid positive coordinate block.')
    if len(proof['pivot_intervals']) != len(prefix)+1 or not all(int(x[0]) > 0 for x in proof['pivot_intervals'][:-1]):
        raise ValueError('Every reused pivot must have a strictly positive lower endpoint.')
    remaining = [n for n in range(N+1) if n not in set(prefix)]
    if not remaining or len(remaining) > args.max_remainder:
        raise ValueError('Unresolved block outside the requested small-matrix scope.')
    bits = packet['grid_bits']
    G = [[interval_ball(x,bits) for x in row] for row in packet['generators']]
    K = arb_mat([[interval_ball(x,bits) for x in row] for row in packet['kernel']])
    kap = interval_ball(packet['kappa'],bits)
    rank = K.nrows()
    if (len(G),len(G[0]),K.ncols()) != (N+1,4+2*rank,rank):
        raise ValueError('Incomplete actual target representation.')
    diagonal = [interval_ball(x,bits) for x in packet['diagonal']]

    # Reconstruct the SAME actual target using its independently supplied diagonal.
    # The far term is the same signed finite moment core, not a PSD assumption.
    target = arb_mat(N+1,N+1)
    for m in range(N+1):
        for n in range(m):
            u,v,t,w = G[m][:4]
            a,b,c,d = G[n][:4]
            value = kap*(u*b-v*a-t*b+w*a-u*d+v*c)/(m*m-n*n)
            target[m,n] = value
            target[n,m] = value
    P = arb_mat([row[4:4+rank] for row in G])
    target -= P*K*P.transpose()
    for n in range(N+1):
        target[n,n] = diagonal[n]
    print('same full target reconstructed; seconds',time.monotonic()-started,flush=True)
    positive = arb_mat([[target[i,j] for j in prefix] for i in prefix])
    coupling = arb_mat([[target[i,j] for j in remaining] for i in prefix])
    small = arb_mat([[target[i,j] for j in remaining] for i in remaining])
    inverse_action = positive.solve(coupling,algorithm='precond')
    small -= coupling.transpose()*inverse_action
    small = (small+small.transpose())/2
    print('validated dense Schur solve completed; remaining',len(remaining),
          'seconds',time.monotonic()-started,flush=True)
    eigenvalues = []
    eigen_error = None
    try:
        spectrum = small.eig(algorithm='rump',multiple=True)
        if len(spectrum) != len(remaining):
            raise ArithmeticError('Incomplete eigenvalue count.')
        eigenvalues = [{'real':envelope(x.real,args.output_bits),
                        'imaginary':envelope(x.imag,args.output_bits)} for x in spectrum]
        if all(x.real > 0 and x.imag.contains(0) for x in spectrum):
            verdict = 'positive-reduced-spectrum'
        elif any(x.real < 0 and x.imag.contains(0) for x in spectrum):
            verdict = 'negative-reduced-spectrum'
        else:
            verdict = 'inconclusive-reduced-spectrum'
    except (ValueError,ArithmeticError) as error:
        eigen_error = str(error)
        verdict = 'inconclusive-reduced-spectrum'
    result = {
        'schema':'ACTUAL_C27_STRUCTURED_DENSE_SCHUR_V1','N':N,'M':M,
        'last_mode':packet['metadata']['last_mode'],'basis':packet['metadata']['basis'],
        'working_bits':args.precision,'grid_bits':args.output_bits,
        'target_packet_sha256':target_hash,'pivot_result_sha256':hashlib.sha256(pivot_raw).hexdigest(),
        'reused_positive_coordinate_count':len(prefix),'remaining_modes':remaining,
        'verdict':verdict,'eigenvalue_intervals':eigenvalues,'eigenvalue_error':eigen_error,
        'reduced_schur_intervals':[[envelope(small[i,j],args.output_bits)
                                    for j in range(small.ncols())] for i in range(small.nrows())],
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_sha256':packet['metadata']['input_sha256'],
        'elapsed_seconds':time.monotonic()-started,
        'python_flint':flint.__version__,'FLINT':flint.__FLINT_VERSION__,
        'conditions':['Same complete actual common-form and exterior-domain premises.',
                      'Same positive actual exterior H and archived W majorizing H^-1.',
                      'Sound reused positive-pivot certificate for this exact complete target.'],
        'scope':'Validated dense Schur/eigenvalue application to the unresolved coordinate block of the same full inverse-residual sufficient target. '
                'A positive reduced spectrum plus the positive prefix would establish this lower target under the stated premises. '
                'A negative reduced spectrum rejects this sufficient target only, not the actual Schur/Weil form. '
                'An inconclusive result certifies neither sign. No target/scalar/moment producer rerun, new general theorem, Lean verification, cofinal positivity or RH result.',
    }
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(verdict,'seconds',time.monotonic()-started,flush=True)
    for value in eigenvalues:
        print('eigenvalue real grid interval',value['real'],flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target',type=Path,required=True)
    parser.add_argument('--pivots',type=Path,required=True)
    parser.add_argument('--precision',type=int,default=768)
    parser.add_argument('--output-bits',type=int,default=640)
    parser.add_argument('--max-remainder',type=int,default=32)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if args.precision < max(640,args.output_bits)+96 or args.out.exists():
        parser.error('Use source/output guard bits and a fresh output path.')
    ctx.prec = args.precision
    run(args)


if __name__ == '__main__':
    main()
