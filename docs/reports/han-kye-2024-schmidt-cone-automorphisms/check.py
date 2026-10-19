#!/usr/bin/env python3
"""Han-Kye (arXiv:2410.13120v1) appendix conjecture: sanity checks of the Lean conventions.

For (m, n, k) in {(3,3,2), (3,4,2), (4,4,2), (4,4,3)}:
  (1) the generators Ad_{s(x)t}, the full transpose and (m = n) the flip send |z><z| with SR z <= k to a
      rank-one projector whose range vector has the same Schmidt rank (200 random z);
  (2) a non-local Ad_V (V a random unitary on C^m (x) C^n) sends some |z><z| with SR z = 1 to a rank-one
      projector of Schmidt rank > k, so it does not preserve S_k;
  (3) Lemma 2-3 route: for a random generator composition Theta, every 2x2 minor of Theta(vec z vec w^T)
      vanishes for random independent z, w (rank-one preservation off the conjugate diagonal).
Exit 0 and RESULT PASS when all hold.
"""
import itertools, sys
import numpy as np
rng = np.random.default_rng(20261019)

def sr(v, m, n): return np.linalg.matrix_rank(v.reshape(m, n), tol=1e-9)
def rand_c(*s): return rng.normal(size=s) + 1j * rng.normal(size=s)
def rand_sr(m, n, r):
    return sum(np.kron(rand_c(m), rand_c(n)) for _ in range(r))

def ad(s, t):
    V = np.kron(s, t); return lambda X: V.conj().T @ X @ V
def transpose(X): return X.T
def flip(m):
    P = np.zeros((m * m, m * m))
    for i in range(m):
        for j in range(m): P[j * m + i, i * m + j] = 1
    return lambda X: P @ X @ P.T

def range_vec(X):
    w, U = np.linalg.eigh((X + X.conj().T) / 2)
    return np.linalg.matrix_rank(X, tol=1e-8), U[:, np.argmax(np.abs(w))]

ok = True
for m, n, k in [(3, 3, 2), (3, 4, 2), (4, 4, 2), (4, 4, 3)]:
    gens = [ad(rand_c(m, m), rand_c(n, n)), transpose] + ([flip(m)] if m == n else [])
    for g in gens:
        for _ in range(200):
            z = rand_sr(m, n, rng.integers(1, k + 1)); Y = g(np.outer(z, z.conj()))
            rk, u = range_vec(Y)
            if rk != 1 or sr(u, m, n) != sr(z, m, n): ok = False; print("FAIL(1)", m, n, k)
    Q, _ = np.linalg.qr(rand_c(m * n, m * n)); adV = lambda X: Q.conj().T @ X @ Q
    worst = max(sr(range_vec(adV(np.outer(z, z.conj())))[1], m, n) for z in [rand_sr(m, n, 1) for _ in range(20)])
    if worst <= k: ok = False; print("FAIL(2)", m, n, k, worst)
    Th = lambda X: gens[0](transpose(gens[-1](X)))
    mx = 0.0
    for _ in range(50):
        z, w = rand_c(m * n), rand_c(m * n); Y = Th(np.outer(z, w))
        mx = max(mx, np.abs(np.linalg.svd(Y, compute_uv=False)[1]))
    if mx > 1e-8: ok = False; print("FAIL(3)", m, n, k, mx)
    print(f"m={m} n={n} k={k} generators_ok worst_SR_under_AdV={worst} max_second_singular={mx:.1e}")
print("RESULT", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
