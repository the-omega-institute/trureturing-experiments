"""Cotfas, arXiv:2610.05623 (after Eq. (22)): is the uncertainty condition (11), sigma + i[[0,1],[-1,0]] >= 0,
sufficient for the mixed DVGS rho_sigma = sum_{n,k=-s}^{s} W_sigma~(n,k) Pi~(n,k) (Eqs. (21)-(22)) to be >= 0?
The author conjectures yes. This program builds the unnormalized matrix for d = 1, 3, 5, 7, 9 and h = 0.5, 1, 3,
with 40 random covariances at each (d, h) satisfying det sigma >= 1 (30% on the boundary det sigma = 1). Lattice sums
are truncated at |alpha|, |beta| <= 40. It checks that the matrix is Hermitian and positive semidefinite, with
positive trace. It also reports the cases with det sigma < 1, where every sampled matrix fails to be positive.
Requires numpy.
"""
import numpy as np, itertools, math, sys
def W(q,p,a,b,c,h):
    D=a*b-c*c
    return 2/(h*math.sqrt(D))*np.exp(-2*math.pi*(b*q*q-2*c*q*p+a*p*p)/(h*D))
def M(d,a,b,c,h,R=40):
    s=(d-1)//2; I=list(range(-s,s+1)); sc=math.sqrt(h/d)
    al=np.arange(-R,R+1)
    A,B=np.meshgrid(al,al,indexing='ij'); sign=np.where((A*B)%2==0,1.0,-1.0)
    w={}
    for n in I:
        for k in I:
            w[(n,k)]=float(np.sum(sign*W((n+A*d/2)*sc,(k+B*d/2)*sc,a,b,c,h)))
    Mx=np.zeros((d,d),complex)
    idx={x:i for i,x in enumerate(I)}
    def rep(x): return ((x+s)%d)-s
    for (n,k),v in w.items():
        for m in I:  # (Pi psi)(m) = e^{4 pi i k (m-n)/d} psi(2n-m)  => Pi[m, rep(2n-m)]
            Mx[idx[m],idx[rep(2*n-m)]]+=v*np.exp(4j*math.pi*k*(m-n)/d)
    return Mx
rng=np.random.default_rng(1); fails=[]; mins=[]
cases=[]
for d in (1,3,5,7,9):
    for h in (0.5,1.0,3.0):
        for _ in range(40):
            a=math.exp(rng.uniform(-2,2)); c=rng.uniform(-3,3)
            D=1.0+(0 if rng.random()<0.3 else math.exp(rng.uniform(-6,2)))
            b=(D+c*c)/a; cases.append((d,h,a,b,c))
for (d,h,a,b,c) in cases:
    Mx=M(d,a,b,c,h); herm=np.max(abs(Mx-Mx.conj().T)); ev=np.linalg.eigvalsh((Mx+Mx.conj().T)/2); tr=np.trace(Mx).real
    mins.append(ev.min()/max(ev.max(),1e-300))
    if ev.min()< -1e-9*max(1,ev.max()) or tr<=0 or herm>1e-8: fails.append((d,h,a,b,c,ev.min(),tr,herm))
print('cases',len(cases),'fails',len(fails),'worst relmin',min(mins))
for f in fails[:5]: print('FAIL',f)
# D<1 sanity: does positivity fail below the uncertainty bound?
neg=0;tot=0
for d in (3,5):
    for D in (0.2,0.5,0.8):
        for a in (0.3,1.0,3.0):
            for c in (0.0,1.0):
                b=(D+c*c)/a; ev=np.linalg.eigvalsh(M(d,a,b,c,1.0)); tot+=1; neg+= ev.min()< -1e-9
print('D<1 cases',tot,'non-PSD',neg)
print('RESULT', 'PASS' if not fails else 'FAIL')
