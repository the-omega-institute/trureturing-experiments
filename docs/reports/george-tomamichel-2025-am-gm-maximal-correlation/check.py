"""Check for the first Section VII question of George and Tomamichel (arXiv:2505.15281v3): mu_AM < 1 and mu_GM = 1?

mu_f is Definition 24 for f_AM(x) = (x+1)/2 (inner product Tr(sigma X^2) on Hermitian X) and f_GM(x) = sqrt x
(Tr(sigma^1/2 X sigma^1/2 X)): the top singular value of the bilinear form Tr[rho (X (x) Y)] on centred Hermitian
matrices in those metrics.  The script checks, numerically: mu_AM <= mu_GM <= 1 on 200 random two-qubit states (the
source states mu_AM <= mu_GM); mu_AM = mu_GM on random classical (diagonal) states (the source's Proposition 31);
mu(rho (x) rho) = mu(rho) for both coefficients on 5 random two-qubit states; and mu_GM = mu_AM = 1 on a state with
perfectly correlated block structure.
"""
import numpy as np
rng = np.random.default_rng(3)
def herm_basis(d):
    B=[]
    for i in range(d):
        E=np.zeros((d,d),complex); E[i,i]=1; B.append(E)
    for i in range(d):
        for j in range(i+1,d):
            E=np.zeros((d,d),complex); E[i,j]=E[j,i]=1/np.sqrt(2); B.append(E)
            F=np.zeros((d,d),complex); F[i,j]=-1j/np.sqrt(2); F[j,i]=1j/np.sqrt(2); B.append(F)
    return B
def sqrtm(A):
    w,V=np.linalg.eigh(A); return V@np.diag(np.sqrt(np.clip(w,0,None)))@V.conj().T
def ptrace(rho,dA,dB,keep):
    r=rho.reshape(dA,dB,dA,dB)
    return np.einsum('ijkj->ik',r) if keep=='A' else np.einsum('ijil->jl',r)
def mu(rho,dA,dB,kind):
    sA=ptrace(rho,dA,dB,'A'); sB=ptrace(rho,dA,dB,'B')
    BA=herm_basis(dA); BB=herm_basis(dB)
    def gram(B,s):
        if kind=='AM': return np.array([[np.trace(s@(X@Y+Y@X)/2).real for Y in B] for X in B])
        r=sqrtm(s); return np.array([[np.trace(r@X@r@Y).real for Y in B] for X in B])
    GA,GB=gram(BA,sA),gram(BB,sB)
    M=np.array([[np.trace(rho@np.kron(X,Y)).real for Y in BB] for X in BA])
    # centered subspaces: Tr(sX)=0
    cA=np.array([np.trace(sA@X).real for X in BA]); cB=np.array([np.trace(sB@Y).real for Y in BB])
    def basis_perp(c,G):
        # orthonormal basis (w.r.t G) of {v : c.v = 0}
        n=len(c); P=np.eye(n)-np.outer(c,c)/ (c@c)
        U,S,_=np.linalg.svd(P); K=U[:, :n-1]
        L=np.linalg.cholesky(K.T@G@K); return K@np.linalg.inv(L).T
    KA=basis_perp(cA,GA); KB=basis_perp(cB,GB)
    return np.linalg.svd(KA.T@M@KB,compute_uv=False)[0]
def rand_state(d):
    G=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)); W=G@G.conj().T; return W/np.trace(W).real
ok = True
worst = 0
for _ in range(200):
    rho = rand_state(4)
    a, g = mu(rho, 2, 2, 'AM'), mu(rho, 2, 2, 'GM')
    worst = max(worst, a - g)
    ok &= g <= 1 + 1e-9
print(f"200 random two-qubit states: max(mu_AM - mu_GM) = {worst:.3e}; all mu_GM <= 1: {ok}")
ok &= worst <= 1e-9
dev = 0
for _ in range(50):
    p = rng.random((3, 3)); p /= p.sum()
    rho = np.diag(p.reshape(-1)).astype(complex)
    dev = max(dev, abs(mu(rho, 3, 3, 'AM') - mu(rho, 3, 3, 'GM')))
print(f"50 random classical 3x3 states: max |mu_AM - mu_GM| = {dev:.3e} (source Proposition 31: equal)")
ok &= dev < 1e-8
for _ in range(5):
    rho = rand_state(4)
    r2 = np.kron(rho, rho).reshape(2, 2, 2, 2, 2, 2, 2, 2).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(16, 16)
    g1, g2 = mu(rho, 2, 2, 'GM'), mu(r2, 4, 4, 'GM')
    a1, a2 = mu(rho, 2, 2, 'AM'), mu(r2, 4, 4, 'AM')
    print(f"mu_GM(rho) = {g1:.6f}, mu_GM(rho x rho) = {g2:.6f}; mu_AM(rho) = {a1:.6f}, mu_AM(rho x rho) = {a2:.6f}")
    ok &= abs(g1 - g2) < 1e-8 and abs(a1 - a2) < 1e-8
a = rand_state(2); b = rand_state(2)
A1 = np.zeros((3, 3), complex); A1[:2, :2] = a; B1 = np.zeros((3, 3), complex); B1[:2, :2] = b
A2 = np.zeros((3, 3), complex); A2[2, 2] = 1; B2 = np.zeros((3, 3), complex); B2[2, 2] = 1
rho = 0.6 * np.kron(A1, B1) + 0.4 * np.kron(A2, B2)
g, a_ = mu(rho, 3, 3, 'GM'), mu(rho, 3, 3, 'AM')
print(f"block state 0.6 a(x)b + 0.4 |2><2|(x)|2><2|: mu_GM = {g:.12f}, mu_AM = {a_:.12f}")
ok &= abs(g - 1) < 1e-9 and abs(a_ - 1) < 1e-9
print("ALL_CHECKS", ok)
raise SystemExit(0 if ok else 1)
