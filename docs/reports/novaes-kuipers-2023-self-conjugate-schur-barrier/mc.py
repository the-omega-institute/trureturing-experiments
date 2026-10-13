import numpy as np, itertools
rng=np.random.default_rng(1)
def haar(M):
    Z=(rng.normal(size=(M,M))+1j*rng.normal(size=(M,M)))/np.sqrt(2)
    Q,R=np.linalg.qr(Z); d=np.diag(R); return Q*(d/abs(d))
def schur(lam,x):
    # bialternant formula
    n=len(x); lam=list(lam)+[0]*(n-len(lam))
    if len(lam)>n: return 0.0
    num=np.array([[xi**(lam[j]+n-1-j) for j in range(n)] for xi in x])
    den=np.array([[xi**(n-1-j) for j in range(n)] for xi in x])
    return np.linalg.det(num)/np.linalg.det(den)
def P(U,r):
    M=U.shape[0]; I=np.eye(M); A=np.linalg.inv(I+r*U)
    return (1-r*r)*A@A.conj().T
for M,lam in [(2,(1,)),(2,(2,1)),(2,(2,)),(3,(2,1)),(3,(2,2)),(3,(3,1,1)),(3,(3,1)),(2,(2,2))]:
    for r in (0.3,0.6):
        N=40000; acc=0
        for _ in range(N):
            ev=np.linalg.eigvalsh(P(haar(M),r)); acc+=schur(lam,ev)
        target=schur(lam,np.ones(M)+1e-7*np.arange(M))
        print(M,lam,r,round(acc/N,3),'target',round(target,3))
