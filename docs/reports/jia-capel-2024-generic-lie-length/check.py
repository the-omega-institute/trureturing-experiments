"""Check of a proposed proof of Jia-Capel, 'A generic quantum Wielandt's inequality' (arXiv:2301.08241v3; Quantum 8
(2024) 1331), Section 6.3, Conjecture 6: a random Lie-generating set S of su(n) has Lie-length Theta(log n) almost surely.

Construction (n >= 2): vertices 0..n-1; binary heap on 1..n-1 (parent(v) = v // 2) plus the anchor edge 0-1;
colours c(1) = 0, c(2v) = c(v) + 1, c(2v+1) = c(v) + 2 (mod 3); H = diag(7, lambda_{c(1)}, ...), lambda = (0, 1, 3);
A = H - tr(H) I / n, B = adjacency matrix; generators iA, iB in su(n).  Claimed: Lie-length <= 16 (1 + floor(log2(n-1))).

(1) Exact combinatorial conditions for Lemma 1 (n = 2..2000): every edge gap H_uu - H_vv lies in
    Delta = {+-1, +-2, +-3, +-7}; gaps +-7 occur only on the anchor edge; siblings and parent-child pairs get distinct
    colours, so for distinct colours (a, b) the gap lambda_b - lambda_a identifies the colour pair.
(2) Exact integer check of the propagation identity of Lemma 2 for n = 2..40: [T_ba, E_v0] = E_w0 + R with R = E_u0
    when the parent u of v has colour b, else 0 (and the transposed identity).
(3) Numerical Lie-length (float rank, tolerance 1e-8) of the construction and of random pairs, n = 2..14:
    lower bound ceil(log2(n^2+1)) - 1 <= L <= 16 (1 + floor(log2(n-1))).
"""
import math
import numpy as np
bad = 0
lam = (0, 1, 3)
def colours(n):
    c = {1: 0}
    for v in range(2, n):
        c[v] = (c[v // 2] + (1 if v % 2 == 0 else 2)) % 3
    return c
def hdiag(n, c): return [7] + [lam[c[v]] for v in range(1, n)]
Delta = {1, 2, 3, 7, -1, -2, -3, -7}
for n in range(2, 2001):
    c = colours(n); h = hdiag(n, c)
    edges = [(0, 1)] + [(v, v // 2) for v in range(2, n)]
    for (u, v) in edges:
        for (x, y) in ((u, v), (v, u)):
            g = h[x] - h[y]
            if g not in Delta: bad += 1
            if abs(g) == 7 and {x, y} != {0, 1}: bad += 1
    for v in range(2, n):
        if c[v] == c[v // 2]: bad += 1
    for v in range(1, n):
        if 2 * v + 1 < n and c[2 * v] == c[2 * v + 1]: bad += 1
gaps = {(a, b): lam[b] - lam[a] for a in range(3) for b in range(3) if a != b}
if len(set(gaps.values())) != 6: bad += 1
print("(1) combinatorial conditions checked n=2..2000, bad so far", bad)

def E(n, i, j):
    M = np.zeros((n, n), dtype=np.int64); M[i, j] = 1; return M
for n in range(3, 41):
    c = colours(n)
    def T(b, a):
        M = np.zeros((n, n), dtype=np.int64)
        for v in range(2, n):
            p = v // 2
            for (x, y) in ((v, p), (p, v)):
                if c[x] == b and c[y] == a: M[x, y] = 1
        return M
    for v in range(1, n):
        for w in (2 * v, 2 * v + 1):
            if w >= n: continue
            a, b = c[v], c[w]
            lhs = T(b, a) @ E(n, v, 0) - E(n, v, 0) @ T(b, a)
            rhs = E(n, w, 0).copy()
            u = v // 2
            if u >= 1 and c[u] == b: rhs += E(n, u, 0)
            if not np.array_equal(lhs, rhs): bad += 1
            lhs2 = E(n, 0, v) @ T(a, b) - T(a, b) @ E(n, 0, v)
            rhs2 = E(n, 0, w).copy()
            if u >= 1 and c[u] == b: rhs2 += E(n, 0, u)
            if not np.array_equal(lhs2, rhs2): bad += 1
print("(2) propagation identities checked n=3..40, bad so far", bad)

rng = np.random.default_rng(6)
def lie_length(S, n, maxL=60):
    target = n * n - 1; basis = []
    def add(Ms):
        for M in Ms:
            v = np.concatenate([M.real.ravel(), M.imag.ravel()])
            for b in basis: v = v - np.dot(b, v) * b
            nv = np.linalg.norm(v)
            if nv > 1e-8: basis.append(v / nv)
    def indep(Ms):
        out, B = [], []
        for M in Ms:
            v = np.concatenate([M.real.ravel(), M.imag.ravel()])
            for b in B: v = v - np.dot(b, v) * b
            nv = np.linalg.norm(v)
            if nv > 1e-9: B.append(v / nv); out.append(M / np.linalg.norm(M))
        return out
    U = list(S); add(U)
    if len(basis) >= target: return 1
    for L in range(2, maxL + 1):
        U = indep([X @ Y - Y @ X for X in S for Y in U]); add(U)
        if len(basis) >= target: return L
    return None
def construct(n):
    c = colours(n); H = np.diag(hdiag(n, c)).astype(float); A = H - np.trace(H) / n * np.eye(n)
    B = np.zeros((n, n)); B[0, 1] = B[1, 0] = 1
    for v in range(2, n): B[v, v // 2] = B[v // 2, v] = 1
    return [1j * A, 1j * B]
def rand_su(n):
    X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); X = X - X.conj().T; return X - np.trace(X) / n * np.eye(n)
for n in range(2, 15):
    up = 16 * (1 + math.floor(math.log2(n - 1))); lo = math.ceil(math.log2(n * n + 1)) - 1
    Lc = lie_length(construct(n), n); Lr = lie_length([rand_su(n), rand_su(n)], n)
    if Lc is None or not (lo <= Lc <= up): bad += 1
    if Lr is None or not (lo <= Lr <= up): bad += 1
    print("(3) n", n, "construction", Lc, "random", Lr, "bounds", lo, up)
print("bad=", bad)
