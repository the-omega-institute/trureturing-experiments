"""Numerical checks for Tavakoli-Morelli, arXiv:2402.09972v2 (PRA 110, 062417 (2024)), Conjectures 2 and 3:
for every state of Schmidt number r,
  EAM (Conj. 2): ||P_n||_tr <= (d(d-1) + r(n-d)) / (n(n-1)),  P_ab = Tr[rho (E_a (x) F_b)], E_a = (d/n)|psi_a><psi_a|;
  MUB (Conj. 3): ||Q_m||_tr <= 1 + (m-1) r / d,  Q_(l,a),(k,b) = Tr[rho (E_la (x) F_kb)] with m of the d+1 MUBs.
Random pure states of Schmidt rank r, random local unitaries on both sides' measurements, and mixtures; plus the sharp
MUB family psi_r = r^{-1/2} sum_{i<r} |ii> with Bob's conjugate bases.  Needs numpy.  Prints ALL_OK."""
import numpy as np
rng = np.random.default_rng(11)
ok = True
def mubs(d):  # Wootters-Fields bases for odd prime d: computational + d bases
    w = np.exp(2j*np.pi/d)
    B = [np.eye(d, dtype=complex)]
    for k in range(d):
        M = np.array([[w**((k*j*j + l*j) % d) for l in range(d)] for j in range(d)]) / np.sqrt(d)
        B.append(M)          # columns are the basis vectors
    return B
def haar(d):
    q, r = np.linalg.qr(rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d)))
    return q * (np.diag(r)/abs(np.diag(r)))
def schmidt_state(d, r):
    s = rng.random(r); s /= np.linalg.norm(s)
    psi = np.zeros((d, d), dtype=complex)
    for i in range(r): psi[i, i] = s[i]
    U, V = haar(d), haar(d)
    return (U @ psi @ V.T).reshape(-1)
def rho_mix(d, r, k):
    ps = rng.random(k); ps /= ps.sum(); rho = 0
    for p in ps:
        v = schmidt_state(d, r); rho = rho + p*np.outer(v, v.conj())
    return rho
def table(rho, EA, FB):
    d = EA[0].shape[0]
    return np.array([[np.real(np.trace(rho @ np.kron(E, F))) for F in FB] for E in EA])
def trn(M): return np.linalg.svd(M, compute_uv=False).sum()
# MUB, d = 3, 5
for d in (3, 5):
    B = mubs(d)
    for m in range(1, d+1):
        for r in range(1, d+1):
            worst = 0
            for t in range(60):
                idxA = rng.choice(d+1, m, replace=False); idxB = rng.choice(d+1, m, replace=False)
                UA, UB = haar(d), haar(d)
                EA = [np.outer(UA@B[l][:, a], (UA@B[l][:, a]).conj()) for l in idxA for a in range(d)]
                FB = [np.outer(UB@B[l][:, a], (UB@B[l][:, a]).conj()) for l in idxB for a in range(d)]
                rho = rho_mix(d, r, rng.integers(1, 4))
                worst = max(worst, trn(table(rho, EA, FB)) - (1 + (m-1)*r/d))
            ok &= worst <= 1e-9
            # sharp family: Bob uses the conjugate bases, psi_r with first basis computational
            psi = np.zeros((d, d), dtype=complex); psi[:r, :r] = np.eye(r)/np.sqrt(r); v = psi.reshape(-1)
            rho = np.outer(v, v.conj())
            EA = [np.outer(B[l][:, a], B[l][:, a].conj()) for l in range(m) for a in range(d)]
            FB = [np.outer(B[l][:, a].conj(), B[l][:, a]) for l in range(m) for a in range(d)]
            sharp = trn(table(rho, EA, FB))
            ok &= abs(sharp - (1 + (m-1)*r/d)) < 1e-9
    print(f"MUB d={d}: random states never exceed 1+(m-1)r/d; sharp family attains it for all m, r: {ok}")
# EAMs: trine (d=2,n=3), tetrahedral SIC (d=2,n=4), real simplex ETF (d=3,n=4)
def eam(vecs, d):
    n = len(vecs); E = [d/n*np.outer(v, v.conj()) for v in vecs]
    S = sum(E); assert np.allclose(S, np.eye(d))
    for i in range(n):
        for j in range(n):
            if i != j: assert abs(abs(np.vdot(vecs[i], vecs[j]))**2 - (n-d)/(d*(n-1))) < 1e-9
    return E
trine = [np.array([np.cos(2*np.pi*k/3), np.sin(2*np.pi*k/3)], dtype=complex) for k in range(3)]
tet = []
for b in ([1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]):
    bx, by, bz = np.array(b)/np.sqrt(3)
    th = np.arccos(bz); ph = np.arctan2(by, bx)
    tet.append(np.array([np.cos(th/2), np.exp(1j*ph)*np.sin(th/2)]))
simplex = [np.array(v, dtype=complex)/np.sqrt(3) for v in ([1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1])]
for name, vecs, d in (("trine", trine, 2), ("tetra-SIC", tet, 2), ("simplex-ETF", simplex, 3)):
    n = len(vecs); E0 = eam(vecs, d)
    for r in range(1, d+1):
        worst = 0
        for t in range(400):
            UA, UB = haar(d), haar(d)
            EA = [UA@E@UA.conj().T for E in E0]; FB = [UB@E@UB.conj().T for E in E0]
            rho = rho_mix(d, r, rng.integers(1, 4))
            worst = max(worst, trn(table(rho, EA, FB)) - (d*(d-1) + r*(n-d))/(n*(n-1)))
        print(f"EAM {name} d={d} n={n} r={r}: max excess over Conj.2 bound = {worst:.3e}")
        ok &= worst <= 1e-9
print("ALL_OK" if ok else "FAIL")
