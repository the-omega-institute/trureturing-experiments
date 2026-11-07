#!/usr/bin/env python3
"""Numerical diagnostics only: exact-law Gaussian generator, separate from check.py."""
import json
import numpy as np

rng = np.random.default_rng(261039674)
N = 100000
rows = []
for D, M, d in [(1,M,1) for M in range(1,5)] + [(2,1,1)] + [(2,M,2) for M in range(2,5)] + [(20,30,20)]:
    trials = N if D <= 2 else 5000
    p=M*(D+1)
    B=np.zeros((p,p))
    for i in range(M):
        B[i*D,M*D+i]=B[M*D+i,i*D]=0.5
    strict=0; negative=0; misses=0
    c=12
    for start in range(0,trials,250):
        A=rng.standard_normal((250,p,d))
        Q=np.linalg.qr(A,mode='reduced')[0]
        g=rng.standard_normal((250,p))
        Qw=Q[:,np.arange(M)*D,:]; Qa=Q[:,M*D+np.arange(M),:]
        C=np.einsum('nmi,nmj->nij',Qw,Qa)
        R=(C+C.transpose(0,2,1))/2
        ev=np.linalg.eigvalsh(R)
        neg=ev[:,-1] < -1e-10
        small=(g*g).sum(axis=1) < c/3
        strict+=int(((ev[:,-1]<-0.25)&small).sum())
        negative+=int(neg.sum())
        # For negative definite R the exact algebraic maximum is q(g)-b^T R^-1 b.
        # Here it is evaluated numerically, with a strict margin and no certification claim.
        gw=g[:,np.arange(M)*D]; ga=g[:,M*D+np.arange(M)]
        b=(np.einsum('nmi,nm->ni',Qw,ga)+np.einsum('nmi,nm->ni',Qa,gw))/2
        qg=(gw*ga).sum(axis=1)
        inds=np.flatnonzero(neg)
        if len(inds):
            sol=np.linalg.solve(R[inds],b[inds,:,None])[...,0]
            maximum=qg[inds]-(b[inds]*sol).sum(axis=1)
            misses+=int((maximum<c-1e-8).sum())
    assert D > 2 or strict>0, (D,M,d,'no strict sufficient event observed')
    rows.append(dict(D=D,M=M,d=d,p=p,u='12 e_1',samples=trials,
                     negative_restriction=negative,
                     strict_sufficient_miss_event=strict,
                     scalar_supremum_below_c=misses))
print(json.dumps(dict(seed=261039674,role='floating-point diagnostic; not proof',
                     distribution='iid N(0,1) p-by-d columns, Haar span by QR; independent g ~ N(0,I_p)',
                     sufficient_event='lambda_max(Q^T B Q)<-1/4 and ||g||^2<4',rows=rows),indent=2))
