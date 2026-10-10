import numpy as np, sys
from scipy.optimize import minimize
rng=np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 0)
def T(X):
    n=X.shape[0]; I=np.eye(n)
    A=np.kron(I,X)-np.kron(X.T,I)   # vec(XY-YX) = (I⊗X - X^T⊗I) vec(Y)
    return A.conj().T@A
def f(v,n):
    X=(v[:n*n]+1j*v[n*n:]).reshape(n,n); X=X/np.linalg.norm(X)
    w=np.linalg.eigvalsh(T(X))[::-1]
    return -(w[0]+w[2])
n=int(sys.argv[1]); best=0;bx=None
for r in range(60):
    v0=rng.standard_normal(2*n*n)
    res=minimize(f,v0,args=(n,),method='BFGS',options={'maxiter':3000,'gtol':1e-10})
    if -res.fun>best: best=-res.fun; bx=res.x
print(n, 'max lambda1+lambda3 =', repr(best))
X=(bx[:n*n]+1j*bx[n*n:]).reshape(n,n); X=X/np.linalg.norm(X)
w=np.linalg.eigvalsh(T(X))[::-1]; print('top eigs', np.round(w[:6],6)); print('normal defect', np.linalg.norm(X@X.conj().T-X.conj().T@X))
