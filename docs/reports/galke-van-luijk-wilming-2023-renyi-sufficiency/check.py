import numpy as np, itertools
from numpy.linalg import eigh
def fmp(S,g):
    w,V=eigh(S); return V@np.diag(w**g)@V.conj().T
np.set_printoptions(precision=6)
def build(r, s, eps):
    d=2*r+1; A=np.zeros((d,d),complex)
    for j in range(1,r+1):
        a,b=2*j-1,2*j
        A[0,a]=A[a,0]=1; A[a,b]=A[b,a]=1
        A[b,0]=s[j-1]*1j; A[0,b]=np.conj(A[b,0])
    return (np.eye(d)+eps*A)/d
def sandwiched(rho,sig,al):
    g=(1-al)/(2*al); S=fmp(sig,g); M=S@rho@S; w=np.clip(np.linalg.eigvalsh((M+M.conj().T)/2),0,None)
    return np.log((w**al).sum())/(al-1)
def umegaki(rho,sig):
    w,V=eigh(rho); lr=V@np.diag(np.log(w))@V.conj().T; return np.trace(rho@(lr-np.diag(np.log(np.diag(sig))))).real
for r in (2,3):
    d=2*r+1; sig=np.diag(np.arange(1,d+1,dtype=float)); sig/=sig.trace()
    eps=0.9/(2*r)
    pats=list(itertools.product([1,-1],repeat=r))
    base=build(r,pats[0],eps)
    print('r',r,'min eig', min(np.linalg.eigvalsh(build(r,p,eps)).min() for p in pats))
    maxdiff=0
    for p in pats[1:]:
        rp=build(r,p,eps)
        for al in [0.5,0.6,0.75,0.9,1.5,2,3,5,10]:
            maxdiff=max(maxdiff,abs(sandwiched(base,sig,al)-sandwiched(rp,sig,al)))
        maxdiff=max(maxdiff,abs(umegaki(base,sig)-umegaki(rp,sig)))
    print('  max |profile difference| over patterns & alphas:',maxdiff)
    # invariant: product of triangle holonomies tr(rho P0 rho Pa rho Pb)
    def hol(rho,j): return rho[0,2*j-1]*rho[2*j-1,2*j]*rho[2*j,0]
    for p in pats: 
        rp=build(r,p,eps); prod=np.prod([hol(rp,j) for j in range(1,r+1)])
        h=[np.round(hol(rp,j)/abs(hol(rp,j)),3) for j in range(1,r+1)]
        print('  pattern',p,'holonomy phases',h,'pair product',np.round(hol(rp,1)*hol(rp,2)/abs(hol(rp,1)*hol(rp,2)),3))
