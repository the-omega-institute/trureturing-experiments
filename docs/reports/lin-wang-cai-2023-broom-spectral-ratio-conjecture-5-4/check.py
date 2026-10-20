#!/usr/bin/env python3
"""Lin-Wang-Cai, arXiv:2302.10491v1, Conjecture 5.4 (upper bound R_L(T) <= R_L(B_n^{t(n)}), t(n) = (n-3)/2 odd, (n-4)/2 even).

Part A (exact, Fractions only): certify R_L(B_19^7) > R_L(B_19^8) by negative inertia of L - xI at four rational x,
computed by exact symmetric elimination of the full 19x19 matrices, and cross-check with the tridiagonal
leading-minor recurrence of the reduced form.
Part B (floating point): broom scan n = 8..400 (failure set of the broom-restricted statement), exhaustive scan of all
trees n = 8..20 (needs networkx), and the continuum limit of the optimal broom fraction t/n.
usage: check.py [--no-exhaustive]
"""
import sys
from fractions import Fraction as Q

def broom_edges(n, t):
    # Figure 2 labels: center v0, pendants v1..vt, path v_{t+1}..v_{n-1}, edge v0 v_{t+1}
    E = [(0, i) for i in range(1, t + 2)] + [(i, i + 1) for i in range(t + 1, n - 1)]
    assert len(E) == n - 1
    return E

def laplacian(n, E):
    L = [[Q(0)] * n for _ in range(n)]
    for a, b in E:
        L[a][a] += 1; L[b][b] += 1; L[a][b] -= 1; L[b][a] -= 1
    return L

def neg_inertia(M):
    # exact symmetric Gaussian elimination without pivoting; returns (#negative pivots, ok) with ok=False on a zero pivot
    A = [row[:] for row in M]; n = len(A); neg = 0
    for k in range(n):
        p = A[k][k]
        if p == 0:
            return None
        if p < 0:
            neg += 1
        for i in range(k + 1, n):
            f = A[i][k] / p
            if f:
                for j in range(k, n):
                    A[i][j] -= f * A[k][j]
    return neg

def shifted(L, x):
    return [[L[i][j] - (x if i == j else 0) for j in range(len(L))] for i in range(len(L))]

def recurrence_sign_changes(t, a, q=10**6, n=19):
    m = n + 1 - t
    A = [1, q - a, ((t + 1) * q - a) * (q - a) - t * q * q]
    for j in range(3, m + 1):
        d = 2 if j < m else 1
        A.append((d * q - a) * A[-1] - q * q * A[-2])
    assert all(v != 0 for v in A)
    return sum(1 for u, v in zip(A, A[1:]) if (u > 0) != (v > 0))

def part_a():
    n, q = 19, 10**6
    L7, L8 = laplacian(n, broom_edges(n, 7)), laplacian(n, broom_edges(n, 8))
    rows = [("B_19^7", L7, 7, 35467, 2), ("B_19^7", L7, 7, 9018480, 18),
            ("B_19^8", L8, 8, 39388, 1), ("B_19^8", L8, 8, 10014270, 19)]
    ok = True
    for name, L, t, a, want in rows:
        got = neg_inertia(shifted(L, Q(a, q)))
        rec = recurrence_sign_changes(t, a) + (t - 1 if a > q else 0)
        print(f"A {name} x={a}/{q}: negative inertia (full elimination) = {got}, reduced recurrence = {rec}, expected {want}")
        ok &= (got == want == rec)
    # inertia -> eigenvalue bounds (one zero eigenvalue, connected):
    # B^7: #eig < 35467/q = 2 -> mu_{n-1} < 35467/q ; #eig < 9018480/q = 18 -> mu_1 > 9018480/q
    # B^8: #eig < 39388/q = 1 -> mu_{n-1} > 39388/q ; #eig < 10014270/q = 19 -> mu_1 < 10014270/q
    lhs, rhs = 9018480 * 39388, 10014270 * 35467
    print(f"A R_L(B_19^7) > 9018480/35467 and R_L(B_19^8) < 10014270/39388; 9018480*39388 = {lhs} > {rhs} = 10014270*35467: {lhs > rhs}")
    ok &= lhs > rhs
    print("A_OK" if ok else "A_FAIL")
    return ok

def part_b(exhaustive):
    import numpy as np
    from scipy.linalg import eigh_tridiagonal
    from scipy.optimize import brentq, minimize_scalar
    def ratio_tri(n, t):
        m = n - t - 1
        d = [1.0, t + 1.0] + [2.0] * (m - 1) + [1.0]
        e = [-np.sqrt(t)] + [-1.0] * m
        w = np.sort(np.concatenate([eigh_tridiagonal(np.array(d), np.array(e), eigvals_only=True), [1.0] * (t - 1)]))
        return w[-1] / w[1]
    def ratio_full(n, E):
        L = np.zeros((n, n))
        for a, b in E: L[a, a] += 1; L[b, b] += 1; L[a, b] -= 1; L[b, a] -= 1
        w = np.linalg.eigvalsh(L); return w[-1] / w[1]
    conj = lambda n: (n - 3) // 2 if n % 2 else (n - 4) // 2
    assert abs(ratio_tri(19, 7) - ratio_full(19, broom_edges(19, 7))) < 1e-9
    holds, fails = [], []
    for n in range(8, 401):
        rs = {t: ratio_tri(n, t) for t in range(1, n - 1)}
        (fails if max(rs.values()) > rs[conj(n)] * (1 + 1e-12) else holds).append(n)
    print("B broom family, n = 8..400: conjectured broom is the best broom for n in", holds, f"; fails for the other {len(fails)} n")
    for n in (19, 21, 60, 400, 3200):
        ts = range(1, n - 1) if n <= 400 else range(int(0.4 * n), int(0.5 * n))
        r, t = max((ratio_tri(n, t), t) for t in ts)
        print(f"B n={n}: best broom t*={t} (t*/n={t / n:.4f}) ratio {r:.4f}; conjectured t={conj(n)} ratio {ratio_tri(n, conj(n)):.4f}")
    th = lambda mu: brentq(lambda x: np.tan(x) + mu * x, np.pi / 2 + 1e-12, np.pi - 1e-15)
    g = lambda a: a * (1 - a) ** 2 / th(a / (1 - a)) ** 2
    r = minimize_scalar(lambda a: -g(a), bounds=(0.3, 0.6), method="bounded", options={"xatol": 1e-12})
    print(f"B continuum limit: optimal broom fraction alpha* = {r.x:.6f}; asymptotic ratio gain over t = n/2: {g(r.x) / g(0.5):.6f}")
    if exhaustive:
        import networkx as nx
        for n in range(8, 21):
            best = None
            batch, trees = [], []
            def flush():
                nonlocal best, batch, trees
                if batch:
                    w = np.linalg.eigvalsh(np.array(batch)); rr = w[:, -1] / w[:, 1]; i = int(np.argmax(rr))
                    if best is None or rr[i] > best[0]: best = (float(rr[i]), trees[i])
                batch, trees = [], []
            cnt = 0
            for T in nx.nonisomorphic_trees(n):
                E = list(T.edges()); L = np.zeros((n, n))
                for a, b in E: L[a, a] += 1; L[b, b] += 1; L[a, b] -= 1; L[b, a] -= 1
                batch.append(L); trees.append(E); cnt += 1
                if len(batch) >= 20000: flush()
            flush()
            G = nx.Graph(best[1]); iso = [t for t in range(1, n - 1) if nx.is_isomorphic(G, nx.Graph(broom_edges(n, t)))]
            print(f"B all {cnt} trees on n={n}: max R_L = {best[0]:.6f} attained by broom t={iso}; conjectured t={conj(n)}")

if __name__ == "__main__":
    ok = part_a()
    part_b("--no-exhaustive" not in sys.argv)
    print("EXIT=0" if ok else "EXIT=1")
    sys.exit(0 if ok else 1)
