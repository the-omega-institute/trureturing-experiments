"""Fuchs, quant-ph/9601020 (1996 thesis), §3.5.3 after Eq. (3.244): does an orthonormal basis of the Hilbert space make
M(t) = I + t A^{-1/2} B A^{-1/2} diagonal in the induced column-vectorization basis? Here A = I (x) rho0 + conj(rho0) (x) I,
B = I (x) D + conj(D) (x) I, D = rho1 - rho0. Proposed answer: such a basis exists iff [rho0, rho1] = 0.
Numerics (not a proof): minimize the off-diagonal mass of (W^dag M(t) W), W = conj(U) (x) U, over unitaries U (many restarts).
1. Non-commuting qutrit family rho0 = diag(1/2, 1/3, 1/6), rho1 = (1/16)[[7,-1,0],[-1,5,0],[0,0,4]]: minimum stays bounded away
   from zero for t in {0.25, 0.5, 1}.
2. Control: commuting pairs reach (numerically) zero, at the common eigenbasis.
3. Random faithful pairs: the minimum is ~0 iff the commutator vanishes (small commutators give small minima).
"""
import sys
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm, sqrtm

rng = np.random.default_rng(20261017)

def Mt(r0, r1, t):
    d = r0.shape[0]; I = np.eye(d); D = r1 - r0
    A = np.kron(I, r0) + np.kron(r0.conj(), I); B = np.kron(I, D) + np.kron(D.conj(), I)
    w, V = np.linalg.eigh(A); Am = V @ np.diag(w ** -0.5) @ V.conj().T
    return np.eye(d * d) + t * Am @ B @ Am

def unitary(params, d):
    H = np.zeros((d, d), complex); iu = np.triu_indices(d, 1)
    n = len(iu[0])
    H[iu] = params[:n] + 1j * params[n:2 * n]; H = H + H.conj().T
    H[np.diag_indices(d)] = params[2 * n:2 * n + d]
    return expm(1j * H)

def offdiag(M, U):
    W = np.kron(U.conj(), U)
    N = W.conj().T @ M @ W
    return np.linalg.norm(N - np.diag(np.diag(N))) ** 2

def best(M, d, restarts=40):
    k = d * (d - 1) + d; b = np.inf
    for _ in range(restarts):
        res = minimize(lambda p: offdiag(M, unitary(p, d)), rng.normal(size=k) * 2, method='BFGS', options={'gtol': 1e-12, 'maxiter': 4000})
        b = min(b, res.fun)
    return b

def rand_state(d):
    G = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)); r = G @ G.conj().T + 0.2 * np.eye(d); return r / np.trace(r).real

def main():
    fails = []
    r0 = np.diag([1 / 2, 1 / 3, 1 / 6]).astype(complex); r1 = np.array([[7, -1, 0], [-1, 5, 0], [0, 0, 4]], complex) / 16
    comm = np.linalg.norm(r0 @ r1 - r1 @ r0)
    for t in (0.25, 0.5, 1.0):
        M = Mt(r0, r1, t); ev = np.linalg.eigvalsh((M + M.conj().T) / 2)
        b = best(M, 3)
        print(f"non-commuting family, t={t}: min eig M(t) = {ev.min():.4f}; min off-diagonal mass over U(3) = {b:.3e} (commutator norm {comm:.3e})")
        if b < 1e-6: fails.append(("family", t))
    # control: commuting pair (common eigenbasis rotated by a random unitary)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
    c0 = Q @ np.diag([0.5, 0.3, 0.2]) @ Q.conj().T; c1 = Q @ np.diag([0.2, 0.45, 0.35]) @ Q.conj().T
    bc = best(Mt(c0, c1, 0.7), 3)
    print(f"control commuting pair: min off-diagonal mass = {bc:.3e} (expected ~0)")
    if bc > 1e-8: fails.append("control")
    rows = []
    for _ in range(12):
        a = rand_state(3); b2 = rand_state(3)
        s = rng.choice([0.0, 1e-3, 1e-1, 1.0])
        X = rand_state(3)
        # make b2 commute with a up to a perturbation of size s
        w, V = np.linalg.eigh(a); b2 = V @ np.diag(np.diag(V.conj().T @ b2 @ V).real) @ V.conj().T + s * (X - np.trace(X) / 3 * np.eye(3))
        b2 = (b2 + b2.conj().T) / 2
        if np.linalg.eigvalsh(b2).min() <= 0: continue
        b2 = b2 / np.trace(b2).real
        cn = np.linalg.norm(a @ b2 - b2 @ a); m = best(Mt(a, b2, 0.5), 3, restarts=15)
        rows.append((cn, m))
    for cn, m in sorted(rows):
        print(f"   random pair: ||[rho0, rho1]|| = {cn:.2e} -> min off-diagonal mass {m:.2e}")
        if cn < 1e-12 and m > 1e-8: fails.append("random-commuting")
        if cn > 1e-2 and m < 1e-8: fails.append("random-noncommuting")
    print("RESULT", "PASS" if not fails else f"FAIL {fails}")
    sys.exit(0 if not fails else 1)

if __name__ == "__main__":
    main()
