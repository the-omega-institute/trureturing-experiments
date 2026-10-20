import numpy as np, itertools
def sqrtm(M):
    w,V=np.linalg.eigh((M+M.conj().T)/2); w=np.clip(w,0,None); return (V*np.sqrt(w))@V.conj().T
def swap(D):
    F=np.zeros((D*D,D*D))
    for i in range(D):
        for j in range(D): F[j*D+i, i*D+j]=1
    return F
def W(X,D):
    S=(np.eye(D*D)+swap(D))/2
    return 2*S@np.kron(X,np.eye(D))@S/(D+1)
def ptrace_second(M,d1,d2):  # trace out second tensor factor
    return np.einsum('ijkj->ik', M.reshape(d1,d2,d1,d2))
def ptrace_first(M,d1,d2):
    return np.einsum('ijil->jl', M.reshape(d1,d2,d1,d2))
def fid(a,b):
    s=sqrtm(a); return np.real(np.trace(sqrtm(s@b@s)))**2
d=r=2; D=d*r
def P(rho): return W(np.kron(rho,np.eye(r)/r),D)          # on (A1E1)(A2E2)
def T(rho):
    Wd=W(rho,d); B=2*(np.eye(r*r)+swap(r))/2/(r*(r+1))
    M=np.kron(Wd,B)  # order A1 A2 E1 E2
    # reorder to A1 E1 A2 E2
    M=M.reshape(d,d,r,r,d,d,r,r).transpose(0,2,1,3,4,6,5,7).reshape(D*D,D*D)
    return M
def sym_ok(M,D):
    F=swap(D); return np.allclose(F@M,M) and np.allclose(M@F,M)
rng=np.random.default_rng(1)
def rand_rho():
    G=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)); R=G@G.conj().T; return R/np.trace(R)
for name,ch in [('P',P),('T',T)]:
    worst=0
    for k in range(400):
        rho=rand_rho() if k else np.diag([1,0]).astype(complex)
        out=ch(rho); assert sym_ok(out,D), name
        assert np.isclose(np.trace(out),1)
        m1=ptrace_second(out,D,D); m2=ptrace_first(out,D,D)
        g=np.kron(rho,np.eye(r)/r)
        inf=1-min(fid(g,m1),fid(g,m2))
        worst=max(worst,inf)
    print(name,'worst one-site infidelity with gamma=rho(x)I/r over samples',worst)
# exact instance rho=diag(3/4,1/4)
rho=np.diag([0.75,0.25]).astype(complex)
print('P marginal',np.round(np.real(np.diag(ptrace_second(P(rho),D,D))),6), 'T marginal', np.round(np.real(np.diag(ptrace_second(T(rho),D,D))),6))
g=np.kron(rho,np.eye(2)/2)
print('fidelities', fid(g,ptrace_second(P(rho),D,D)), (23+np.sqrt(273))/40, fid(g,ptrace_second(T(rho),D,D)), (7+2*np.sqrt(6))/12)
# pure input lower bound check for P: rho=|0><0|, gamma=|0><0|(x)sigma random sigma -> fidelity <= 0.8
rho=np.diag([1,0]).astype(complex); m=ptrace_second(P(rho),D,D); best=0
for k in range(3000):
    s=rand_rho(); g=np.kron(rho,s); best=max(best,fid(g,m))
print('P pure input: max fidelity over gamma=|0><0|(x)sigma', best, '(bound 0.8)')
