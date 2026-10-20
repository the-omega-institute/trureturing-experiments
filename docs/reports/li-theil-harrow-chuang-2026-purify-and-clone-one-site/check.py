"""Li-Theil-Harrow-Chuang (arXiv:2605.21457v2) S3 A conjecture, n = 1, m = 2, one-site worst-case squared-fidelity risk.
P = W_{dr}(rho (x) I_r/r) (purify-and-clone); T = R[W_d(rho) (x) B_r]R^dagger with B_r = 2 S_r / (r(r+1)) (competitor,
CPTP, output supported on Sym^2 of the two AE pairs). Both one-site marginals have the form tau (x) I_r/r, and
max over gamma in conv Pur_r(rho) of F(gamma, tau (x) I_r/r) = F(rho, tau) (monotonicity + Weyl twirl). So the risks are
max over rank<=r rho of 1 - F(rho, tau(rho)), with tau the A-marginal. Checks, for (d,r) in {(2,2),(3,2),(3,3),(4,2)}:
  1. the A-marginals of P and T equal the depolarized eta_{dr} and eta_d forms (exact matrix identity on random inputs);
  2. sampled worst-case risk over random rank<=r states and pure states matches r(d-1)/[2(dr+1)] for P and
     (d-1)/[2(d+1)] for T, and T < P strictly when r >= 2;
  3. T's output is supported on the totally symmetric subspace: S X S = X with S = (I + F_{dr})/2 after reordering.
"""
import numpy as np
rng=np.random.default_rng(20261020)
def swap(D):
    F=np.zeros((D*D,D*D))
    for i in range(D):
        for j in range(D): F[j*D+i, i*D+j]=1
    return F
def W(X,D):
    S=(np.eye(D*D)+swap(D))/2; return 2*S@np.kron(X,np.eye(D))@S/(D+1)
def sqrtm(M):
    w,V=np.linalg.eigh((M+M.conj().T)/2); w=np.clip(w,0,None); return (V*np.sqrt(w))@V.conj().T
def fid(a,b):
    s=sqrtm(a); return np.real(np.trace(sqrtm(s@b@s)))**2
def ptr2(M,d1,d2): return np.einsum('ijkj->ik', M.reshape(d1,d2,d1,d2))
def reorder(d,r):
    # A1 A2 E1 E2 -> A1 E1 A2 E2 permutation matrix on (d*d*r*r)
    n=d*d*r*r; P=np.zeros((n,n))
    for a1 in range(d):
        for a2 in range(d):
            for e1 in range(r):
                for e2 in range(r):
                    src=((a1*d+a2)*r+e1)*r+e2; dst=((a1*r+e1)*d+a2)*r+e2; P[dst,src]=1
    return P
def rand_state(d,k):
    G=rng.normal(size=(d,k))+1j*rng.normal(size=(d,k)); R=G@G.conj().T; return R/np.trace(R)
ok=True
for d,r in ((2,2),(3,2),(3,3),(4,2)):
    D=d*r; Sr=(np.eye(r*r)+swap(r))/2; Br=2*Sr/(r*(r+1)); Rm=reorder(d,r)
    def Pch(rho): return W(np.kron(rho,np.eye(r)/r),D)
    def Tch(rho): return Rm@np.kron(W(rho,d),Br)@Rm.T
    etaD=(D+2)/(2*(D+1)); etad=(d+2)/(2*(d+1))
    worstP=worstT=0.0; marg_ok=sym_ok=True
    SD=(np.eye(D*D)+swap(D))/2
    for trial in range(60):
        k=1 if trial%3==0 else rng.integers(1,r+1); rho=rand_state(d,k)
        mP=ptr2(Pch(rho),D,D); mT=ptr2(Tch(rho),D,D)
        tauP=etaD*rho+(1-etaD)*np.eye(d)/d; tauT=etad*rho+(1-etad)*np.eye(d)/d
        marg_ok &= np.allclose(mP,np.kron(tauP,np.eye(r)/r)) and np.allclose(mT,np.kron(tauT,np.eye(r)/r))
        X=Tch(rho); sym_ok &= np.allclose(SD@X@SD, X)
        worstP=max(worstP,1-fid(rho,tauP)); worstT=max(worstT,1-fid(rho,tauT))
    fP=r*(d-1)/(2*(d*r+1)); fT=(d-1)/(2*(d+1))
    good = marg_ok and sym_ok and abs(worstP-fP)<1e-7 and abs(worstT-fT)<1e-7 and fT<fP
    ok &= good
    print(f"(d,r)=({d},{r}): marginals_ok={marg_ok} T_sym_support={sym_ok} worstP={worstP:.10f} formula={fP:.10f} worstT={worstT:.10f} formula={fT:.10f} gap={fP-fT:.6f} {'ok' if good else 'FAIL'}")
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
