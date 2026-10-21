"""Hutton-Mehta-Vukovic, arXiv:2508.10664v1, Conjecture 1.5. For density matrices sigma_1..sigma_n in M_d(C) and orthonormal
u_1..u_k in C^n (2 <= k <= n), A := sum_{r != s} Tr(tau_r tau_s) >= min_{|S| = k} sum_{i != j in S} Tr(sigma_i sigma_j),
tau_r = sum_i |<i|u_r>|^2 sigma_i. Checks the identities of the proof (Cauchy-Binet weights p_S = |det U_S|^2 sum to 1, and
A - sum_S p_S f(S) = sum_{r != s} ||T_rs||_HS^2) on random instances, and searches for violations by local minimization."""
import itertools
import numpy as np
rng = np.random.default_rng(20250814)

def rand_state(d):
    G = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)); r = G @ G.conj().T; return r / np.trace(r).real

def rand_U(n, k):
    G = rng.normal(size=(n, k)) + 1j * rng.normal(size=(n, k)); Q, _ = np.linalg.qr(G); return Q

def lhs(sig, U):
    n, k = U.shape; W = np.abs(U) ** 2
    tau = [sum(W[i, r] * sig[i] for i in range(n)) for r in range(k)]
    return sum(np.trace(tau[r] @ tau[s]).real for r in range(k) for s in range(k) if r != s)

def rhs(sig, k):
    n = len(sig); G = np.array([[np.trace(a @ b).real for b in sig] for a in sig])
    return min(sum(G[i, j] for i in S for j in S if i != j) for S in itertools.combinations(range(n), k)), G

worst = np.inf; ident_err = 0.0; psum_err = 0.0; trials = 0
for n in range(2, 7):
    for k in range(2, n + 1):
        for d in (1, 2, 3):
            for t in range(40):
                sig = [rand_state(d) for _ in range(n)]
                U = rand_U(n, k)
                A = lhs(sig, U); m, G = rhs(sig, k)
                worst = min(worst, A - m); trials += 1
                pS = {S: abs(np.linalg.det(U[list(S), :])) ** 2 for S in itertools.combinations(range(n), k)}
                psum_err = max(psum_err, abs(sum(pS.values()) - 1))
                Ef = sum(p * sum(G[i, j] for i in S for j in S if i != j) for S, p in pS.items())
                T = {(r, s): sum(np.conj(U[i, r]) * U[i, s] * sig[i] for i in range(n)) for r in range(k) for s in range(k) if r != s}
                hs = sum(np.trace(M.conj().T @ M).real for M in T.values())
                ident_err = max(ident_err, abs((A - Ef) - hs))
print(f"random trials={trials} min(A - min_S f(S))={worst:.3e} max|sum p_S - 1|={psum_err:.2e} max|A - E f - sum||T||^2|={ident_err:.2e}")

# Local minimization of A - min f over U (parametrized by QR of a free matrix) for adversarial states.
def search(n, k, d, iters=300):
    sig = [rand_state(d) for _ in range(n)]
    m, _ = rhs(sig, k)
    X = rng.normal(size=(n, k)) + 1j * rng.normal(size=(n, k)); best = np.inf; step = 0.3
    def f(X):
        Q, _ = np.linalg.qr(X); return lhs(sig, Q) - m
    cur = f(X)
    for _ in range(iters):
        Y = X + step * (rng.normal(size=X.shape) + 1j * rng.normal(size=X.shape))
        v = f(Y)
        if v < cur: X, cur = Y, v
        else: step *= 0.98
    return cur
mins = [search(n, k, d) for n in (3, 4, 5) for k in range(2, n + 1) for d in (1, 2) for _ in range(5)]
print(f"local search runs={len(mins)} min(A - min f)={min(mins):.3e}")
