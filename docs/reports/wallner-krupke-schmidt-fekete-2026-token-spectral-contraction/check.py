# Wallner–Krupke–Schmidt–Fekete arXiv:2608.27085 Conjecture 1: rho(A^{-1} B_alpha) < 1 (Token-1).
import numpy as np, random
def mats(v, alpha):
    n=len(v); m=n-1
    A=np.zeros((m,m)); B=np.zeros((m,m))
    for i in range(1,m+1):          # 1-based pair index, robots i, i+1 (v[i-1], v[i])
        r=i-1; vi=v[i-1]; vj=v[i]
        if i%2==1:
            A[r,r]=1.0
            B[r,r]=-1.0
            if i+1<=m: B[r,r+1]=2*vi/(vi+vj)
            if i-1>=1: B[r,r-1]=2*vj/(vi+vj)
        else:
            A[r,r]=vi+vj
            if i-1>=1: A[r,r-1]=-2*vj
            if i+1<=m: A[r,r+1]=-2*vi
            B[r,r]=-(vi+vj)
    Ba=B.copy()
    v1,v2=v[0],v[1]
    Ba[0,:]=0.0
    Ba[0,0]=(v2*(1-2*alpha)-alpha*v1)/(v2+alpha*v1)
    if m>=2: Ba[0,1]=2*alpha*v1/(v2+alpha*v1)
    return A,B,Ba
random.seed(7); worst=0; tests=0; undamped_dev=0
for n in range(2,13):
    for t in range(400):
        v=[random.uniform(0.05,20) if random.random()<0.5 else random.choice([1,2,3,0.5]) for _ in range(n)]
        alpha=random.choice([random.uniform(1e-3,0.999),0.5,0.01,0.99])
        A,B,Ba=mats(v,alpha)
        M=np.linalg.solve(A,B); Ma=np.linalg.solve(A,Ba)
        undamped_dev=max(undamped_dev,np.max(np.abs(np.abs(np.linalg.eigvals(M))-1)))
        r=np.max(np.abs(np.linalg.eigvals(Ma))); worst=max(worst,r); tests+=1
        # paper's alpha=1 consistency: Ba == B
        _,_,B1=mats(v,1.0); assert np.allclose(B1,B)
print("tests",tests,"max rho(M_alpha)",worst,"max | |eig(M)| - 1 |",undamped_dev)
