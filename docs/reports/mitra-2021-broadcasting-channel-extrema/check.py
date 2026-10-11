#!/usr/bin/env python3
"""Checks for A. Mitra, "Layers of classicality in the compatibility of measurements",
arXiv:2101.05752v4 (Phys. Rev. A 104, 022206), Section III.A: does the set of channels broadcasting
an observable A have a highest / a lowest element under nondisturbing precomposition?
Usage: python3 check.py (numpy).  Exits 0 and prints `bad= 0` when every check passes."""
import sys
import numpy as np

bad = 0
def ok(cond, msg):
    global bad
    if not cond:
        bad += 1
        print("FAIL", msg)

rng = np.random.default_rng(20261011)
TOL = 1e-10

def family(r, alpha=0.75, beta=0.25):
    """dim H = r+1, P = first r basis vectors, Q = last; A0 = alpha P + beta Q, A1 = I - A0.
    Output vectors: e_j (x) e_j (j < r) and (e0 e1 + e1 e0)/sqrt2; w = e_r (x) e_r."""
    d = r + 1
    I = np.eye(d)
    P = np.diag([1.0] * r + [0.0]); Q = I - P
    A = [alpha * P + beta * Q, I - (alpha * P + beta * Q)]
    def ket(i, j):
        v = np.zeros(d * d); v[i * d + j] = 1; return v
    vs = [ket(j, j) for j in range(r)] + [(ket(0, 1) + ket(1, 0)) / np.sqrt(2)]
    w = ket(r, r)
    return d, P, Q, A, vs, w

def kraus_of(d, v, w, r):
    """Kraus operators of Lambda_v(X) = Tr(PX) vv* + Tr(QX) ww*: v e_j^* (j<r) and w e_r^*."""
    Ks = [np.outer(v, np.eye(d)[j]) for j in range(r)] + [np.outer(w, np.eye(d)[r])]
    return Ks

def apply(Ks, X):
    return sum(K @ X @ K.conj().T for K in Ks)

def ptrace(M, d, keep):
    T = M.reshape(d, d, d, d)
    return np.einsum('ijkj->ik', T) if keep == 0 else np.einsum('jijk->ik', T)

def rand_state(d):
    z = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    rho = z @ z.conj().T
    return rho / np.trace(rho).real

def rand_unitary(n):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    q, rr = np.linalg.qr(z)
    return q * (np.diag(rr) / abs(np.diag(rr)))

for r in range(2, 6):
    d, P, Q, A, vs, w = family(r)
    # outputs: projectors, mutually orthogonal, orthogonal to w
    S = [np.outer(v, v.conj()) for v in vs]
    for i, s in enumerate(S):
        ok(np.allclose(s @ s, s), f"r={r} sigma_{i} projector")
        ok(abs(np.vdot(w, vs[i])) < TOL, f"r={r} sigma_{i} orthogonal to w")
        for j in range(i):
            ok(abs(np.vdot(vs[i], vs[j])) < TOL, f"r={r} sigma orthogonal {i},{j}")
    ok(len(S) == r + 1 and np.trace(P) == r, f"r={r} counts r+1 > r")
    for i, v in enumerate(vs):
        Ks = kraus_of(d, v, w, r)
        ok(np.allclose(sum(K.conj().T @ K for K in Ks), np.eye(d)), f"r={r} Kraus completeness {i}")
        for _ in range(20):
            rho = rand_state(d)
            out = apply(Ks, rho)
            for x in range(2):
                p = np.trace(rho @ A[x]).real
                ok(abs(np.trace(ptrace(out, d, 0) @ A[x]).real - p) < 1e-9, f"r={r} marginal A {i}")
                ok(abs(np.trace(ptrace(out, d, 1) @ A[x]).real - p) < 1e-9, f"r={r} marginal B {i}")
        # Lemma 1: Lambda o Theta = Lambda for A-nondisturbing Theta.  Random nondisturbing
        # channels: block unitaries on ran P and ran Q mixed with a block measure-and-prepare.
        for _ in range(10):
            Up = rand_unitary(r); U = np.zeros((d, d), complex); U[:r, :r] = Up; U[r, r] = np.exp(1j * rng.normal())
            sig_p = rand_state(r); sigP = np.zeros((d, d), complex); sigP[:r, :r] = sig_p
            sigQ = np.zeros((d, d), complex); sigQ[r, r] = 1
            lam = rng.uniform()
            def theta(X):
                return lam * U @ X @ U.conj().T + (1 - lam) * (np.trace(P @ X) * sigP + np.trace(Q @ X) * sigQ)
            # it is nondisturbing: dual fixes P and Q, hence A
            for _ in range(5):
                rho = rand_state(d)
                for x in range(2):
                    ok(abs(np.trace(theta(rho) @ A[x]) - np.trace(rho @ A[x])) < 1e-9, f"r={r} theta nondisturbing")
                ok(np.allclose(apply(Ks, theta(rho)), apply(Ks, rho), atol=1e-9), f"r={r} Lemma 1 {i}")
    # distinct channels: they differ on e0 e0*
    E00 = np.zeros((d, d)); E00[0, 0] = 1
    outs = [apply(kraus_of(d, v, w, r), E00) for v in vs]
    for i in range(len(outs)):
        for j in range(i):
            ok(not np.allclose(outs[i], outs[j]), f"r={r} channels differ {i},{j}")
    # Lemma 2 bookkeeping: tau supported on ran P has tau <= P
    for _ in range(20):
        t = rand_state(r); tau = np.zeros((d, d), complex); tau[:r, :r] = t
        ok(np.linalg.eigvalsh(P - tau).min() > -1e-12, f"r={r} tau <= P")

# the r = 2 instance of the preregistration with alpha = 3/4, beta = 1/4 also for other alpha != beta
for alpha, beta in ((0.6, 0.3), (0.2, 0.9), (0.75, 0.25)):
    d, P, Q, A, vs, w = family(2, alpha, beta)
    ok(np.allclose(A[0], (beta) * np.eye(d) + (alpha - beta) * P), "A0 = beta I + (alpha-beta) P")

print("bad=", bad)
sys.exit(1 if bad else 0)
