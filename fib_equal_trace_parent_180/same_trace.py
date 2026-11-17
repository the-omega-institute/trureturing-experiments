#!/usr/bin/env python3
"""Exact exceptional-stratum elimination for MM pairs with equal first trace.

This script imposes no condition on the independently selected MJ pair.
It is a partial-stratum calculation, not a dimension certificate for all Bad.
"""
import json
import sys
import time
from pathlib import Path

import flint

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
from incidence import model, encode, info, dump


def main():
    c=flint.fmpq_mpoly_ctx.get(('u','b','e','j','theta','xi'))
    u,b,e,j,theta,xi=c.gens()
    m=model(c,theta,xi)
    dump('universal_0_input.json',{'modulus':0,
        'measurements':{k:info(m[k]) for k in ('R','E','J','Nd','Dd','Ad','phi')},
        'polynomials':{k:encode(m[k]) for k in ('R','E','J','Nd','Dd','Ad','phi')}})
    K=205*u**4+313*u**3-1211*u**2-894*u+1940
    h=u*(u+1)
    M=[[3*h-10,8*u*u+18*u-30,-23*h],
       [10*u*u+8*u-32,10*u*u+19*u-36,-29*h],
       [135*u*u+93*u-420,160*u*u+272*u-560,-381*h]]
    det=(M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
        -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
        +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
    assert det==-20*h*K
    result={'quartic':encode(K),'determinant_identity':True,
            'inputs':{x:info(m[x]) for x in ('Nd','Dd')}}
    dump('same_trace_input.json',result)
    print(json.dumps(result['inputs']),flush=True)
    t=time.monotonic()
    S=m['Nd'].resultant(m['Dd'],1)
    print(json.dumps({'stage':'same_trace_b_resultant','seconds':time.monotonic()-t,'info':info(S)}),flush=True)
    co,fs=S.factor()
    result['b_resultant_seconds']=time.monotonic()-t
    result['b_resultant']=encode(S)
    result['factors']=[{'polynomial':encode(f),'multiplicity':n,'info':info(f)} for f,n in fs]
    dump('same_trace_elimination.json',result)
    print(json.dumps({'factors':[(str(f) if len(f.to_dict())<12 else info(f),n) for f,n in fs]}),flush=True)
    # A nonzero norm of the entire b resultant gives a global polynomial
    # necessary condition on the equal-u, A != 0 stratum.
    t=time.monotonic()
    Z=K.resultant(S,0)
    assert not Z.is_zero()
    co,fs=Z.factor()
    result['parent_resultant_seconds']=time.monotonic()-t
    result['parent_polynomial']=encode(Z)
    result['parent_factors']=[{'polynomial':encode(f),'multiplicity':n,'info':info(f)} for f,n in fs]
    dump('same_trace_elimination.json',result)
    print(json.dumps({'stage':'same_trace_parent_resultant','seconds':time.monotonic()-t,'info':info(Z),
       'factors':[(str(f) if len(f.to_dict())<12 else info(f),n) for f,n in fs]}),flush=True)


if __name__=='__main__':
    main()
