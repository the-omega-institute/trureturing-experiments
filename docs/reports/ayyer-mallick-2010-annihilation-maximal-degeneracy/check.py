#!/usr/bin/env python3
"""Check of the proposed refutation of the Ayyer-Mallick maximal-degeneracy conjecture
(Ayyer-Mallick, arXiv:0910.0693v2, J. Phys. A 43 (2010) 045003, after Eq. (79): "each eigenvalue seems
to have exactly one eigenvector"; Ayyer-Strehl, FPSAC 2010, arXiv:1008.2134, Conjecture 1).

Model (lambda = 1): bits b_1..b_L; bulk 10 -> 01 rate 1, 11 -> 00 rate 1; site 1 flips 0 <-> 1 at
rate alpha; site L: 1 -> 0 at rate beta. M[new, old] = rate; columns sum to zero.

 1. Exact (rational): the functions P_r (1 <= r <= floor((L-1)/2)) satisfy (M^T + 2) P_r = 0 for
    L = 5, 6, 7 at several (alpha, beta) with alpha != 1/2, and are linearly independent (exact rank).
    Hence dim ker(M + 2) = dim ker(M^T + 2) >= floor((L-1)/2) >= 2.
 2. Exact: -2 differs from every other eigenvalue of Theorem 4 at these parameters.
 3. Float (numpy SVD, corroboration only): geometric multiplicity of every eigenvalue, L = 4..8,
    including the alpha = 1/2 control.
"""
import itertools, sys
from fractions import Fraction as F
import numpy as np

def generator(L, a, b, num=F):
    states = list(itertools.product([0, 1], repeat=L)); idx = {s: k for k, s in enumerate(states)}
    N = len(states); M = [[num(0)] * N for _ in range(N)]
    def add(o, n, r):
        M[idx[n]][idx[o]] += r; M[idx[o]][idx[o]] -= r
    for s in states:
        for i in range(L - 1):
            if s[i] == 1 and s[i+1] == 0:
                t = list(s); t[i], t[i+1] = 0, 1; add(s, tuple(t), num(1))
            if s[i] == 1 and s[i+1] == 1:
                t = list(s); t[i], t[i+1] = 0, 0; add(s, tuple(t), num(1))
        t = list(s); t[0] = 1 - s[0]; add(s, tuple(t), a)
        if s[L-1] == 1:
            t = list(s); t[L-1] = 0; add(s, tuple(t), b)
    return states, M

def rank(A):
    A = [row[:] for row in A]; n = len(A); m = len(A[0]); r = 0
    for c in range(m):
        p = next((i for i in range(r, n) if A[i][c] != 0), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        for i in range(n):
            if i != r and A[i][c] != 0:
                f = A[i][c] / A[r][c]; A[i] = [x - f*y for x, y in zip(A[i], A[r])]
        r += 1
    return r

def xi(s, i):
    v = 1
    for k in range(i): v *= (1 - 2*s[k])
    return v

def eigenvalues(L, a, b):
    ev = []
    for k in range(L):
        if 2*k <= L-1: ev += [-2*k, -2*k - 2*a - b]
        if 2*k+1 <= L-1: ev += [-(2*k+1) - b, -(2*k+1) - 2*a]
    return ev

ok = True
for (a, b) in [(F(1,3), F(2,5)), (F(2,3), F(7,4)), (F(3,8), F(5,6))]:
    for L in (5, 6, 7):
        states, M = generator(L, a, b); N = len(states)
        d = 1 - 2*a; R = (L - 1) // 2; vecs = []
        for r in range(1, R + 1):
            f = []
            for s in states:
                g = sum((-1)**(j-1) * xi(s, j) * xi(s, 2*r+1-j) for j in range(1, r+1))
                h = sum((-1)**k * d**(2*r-k) * xi(s, 2*r+1-k) for k in range(1, 2*r+1))
                f.append(F(2) * d**(2*r) * g + 2*h + (-1)**r * d**(2*r) - 1)
            Gf = [sum(M[i][j] * f[i] for i in range(N)) for j in range(N)]
            if any(Gf[j] + 2*f[j] != 0 for j in range(N)): ok = False
            vecs.append(f)
        rk = rank(vecs)
        evs = eigenvalues(L, a, b)
        isolated = sum(1 for e in evs if e == -2) == 1
        line = rk == R and R >= 2 and isolated
        ok &= line
        print(f"exact alpha={a} beta={b} L={L}: {R} eigenfunctions of M^T at -2 verified, rank={rk}, "
              f"-2 distinct from the other {2*L-1} eigenvalues={isolated} {'OK' if line else 'FAIL'}")
# float corroboration
for (a, b) in [(1/3, 2/5), (0.5, 2/5)]:
    for L in range(4, 9):
        _, M = generator(L, a, b, num=float); M = np.array(M)
        evs = sorted(set(round(e, 12) for e in eigenvalues(L, a, b)))
        geo = {}
        for mu in evs:
            sv = np.linalg.svd(M - mu*np.eye(len(M)), compute_uv=False)
            geo[mu] = int(np.sum(sv < 1e-8 * sv.max()))
        print(f"float alpha={a} beta={b} L={L}: geometric multiplicity at -2 = {geo.get(-2.0)}, "
              f"total eigenvector dimension = {sum(geo.values())} (2L = {2*L})")
print("ALL_OK" if ok else "SOME_FAIL")
sys.exit(0 if ok else 1)
