"""Run with numpy and scipy installed.
Numerical checks for Chen-Yu, arXiv:1604.05788 (PRA 94, 022307 (2016)), Conjecture 27:
log2 Sch(U) >= K_Ea(U) >= K_E(U), with the two inequalities equalities at the same time.
(1) Lemma 3 (strict assistance at a nonflat maximizing product output) for U = exp(i theta Z(x)Z), small theta.
(2) The interpolation bound H(T M) - H(M) <= 2 ln A_U on random inputs for random 2x2 and 2x3 unitaries,
    A_U = max nuclear norm of product outputs.  Prints ALL_OK."""
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm
rng = np.random.default_rng(7)

def ent(psi, dA):  # von Neumann entropy (nats) across dA : rest, psi a flat vector
    M = psi.reshape(dA, -1); s = np.linalg.svd(M, compute_uv=False)**2; s = s[s > 1e-15]
    return float(-(s*np.log(s)).sum())

def apply_U(U, psi, a, b, ra, rb):  # psi on (A RA) (B RB) ordered A,RA,B,RB; U on A,B
    T = psi.reshape(a, ra, b, rb).transpose(0, 2, 1, 3).reshape(a*b, ra*rb)
    T = U @ T
    return T.reshape(a, b, ra, rb).transpose(0, 2, 1, 3).reshape(-1)

def unit(v): return v/np.linalg.norm(v)

def KE(U, a, b, tries=40):
    ra, rb = a, b
    best = (-1, None)
    def f(z):
        x = unit(z[:2*a*ra].view(complex)); y = unit(z[2*a*ra:].view(complex))
        return -ent(apply_U(U, np.kron(x, y), a, b, ra, rb), a*ra)
    for _ in range(tries):
        z0 = rng.normal(size=2*(a*ra + b*rb))
        r = minimize(f, z0, method='BFGS')
        if -r.fun > best[0]: best = (-r.fun, r.x)
    z = best[1]; x = unit(z[:2*a*ra].view(complex)); y = unit(z[2*a*ra:].view(complex))
    return best[0], np.kron(x, y)

ok = True
# (1) Lemma 3 for exp(i theta ZZ)
Z = np.diag([1, -1]); ZZ = np.kron(Z, Z)
for theta in (0.15, 0.3):
    U = expm(1j*theta*ZZ); a = b = 2
    k, psi0 = KE(U, a, b)
    phi0 = apply_U(U, psi0, a, b, a, b)
    M = phi0.reshape(a*a, b*b); Uu, s, Vh = np.linalg.svd(M); q = s**2; keep = q > 1e-14
    mu = (q[keep]*np.log(q[keep])).sum()
    coeff = np.zeros_like(s); coeff[keep] = s[keep]*(mu - np.log(q[keep]))
    chi = (Uu[:, :len(s)] * coeff) @ Vh[:len(s), :]
    v = float((np.abs(coeff)**2).sum())
    flat = np.ptp(q[keep]) < 1e-9
    gains = []
    for eps in (1e-2, 3e-3, 1e-3):
        phie = unit((M + eps*chi).reshape(-1))          # output vector, ordering A,RA,B,RB
        psie = apply_U(U.conj().T, phie, a, b, a, b)     # pulled-back input U^dagger phi_eps
        gain = ent(phie, a*a) - ent(psie, a*a)
        gains.append(gain - k)
    print(f"theta={theta}: K_E={k/np.log(2):.6f} bits, output nonflat={not flat}, Var={v:.4g}, net gain over K_E (nats) at eps=1e-2,3e-3,1e-3: {gains}")
    ok &= (not flat) and all(g > 0 for g in gains)
# (2) interpolation bound on random inputs
for (a, b) in ((2, 2), (2, 3)):
    U = np.linalg.qr(rng.normal(size=(a*b, a*b)) + 1j*rng.normal(size=(a*b, a*b)))[0]
    ra, rb = a, b
    def nuc(z):
        x = unit(z[:2*a*ra].view(complex)); y = unit(z[2*a*ra:].view(complex))
        out = apply_U(U, np.kron(x, y), a, b, ra, rb).reshape(a*ra, b*rb)
        return -np.linalg.svd(out, compute_uv=False).sum()
    AU = max(-minimize(nuc, rng.normal(size=2*(a*ra + b*rb)), method='BFGS').fun for _ in range(30))
    worst = -1e9
    for _ in range(3000):
        RA, RB = rng.integers(1, 4), rng.integers(1, 4)
        psi = unit(rng.normal(size=a*RA*b*RB) + 1j*rng.normal(size=a*RA*b*RB))
        if rng.random() < 0.5:   # near-product inputs as well
            psi = unit(np.kron(unit(rng.normal(size=a*RA)+1j*rng.normal(size=a*RA)), unit(rng.normal(size=b*RB)+1j*rng.normal(size=b*RB))) + 0.3*psi)
        d = ent(apply_U(U, psi, a, b, RA, RB), a*RA) - ent(psi, a*RA)
        worst = max(worst, d - 2*np.log(AU))
    print(f"{a}x{b} random U: A_U={AU:.6f} (sqrt Sch bound {np.sqrt(min(a*a, b*b)):.4f}); max of [H(UM)-H(M)-2lnA_U] over 3000 inputs = {worst:.4g}")
    ok &= worst <= 1e-6 and AU <= np.sqrt(min(a*a, b*b)) + 1e-6
print("ALL_OK" if ok else "FAIL")
