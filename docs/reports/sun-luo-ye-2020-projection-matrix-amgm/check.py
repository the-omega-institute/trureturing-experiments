#!/usr/bin/env python3
"""Sun-Luo-Ye, Math. Oper. Res. 45 (2020) 233-271 (arXiv:1503.06387v3), Conjecture 1 (Sec. 3.3):
(1/n!) sum_sigma P_sigma(n) ... P_sigma(1)  <=  ((1/n) sum_i P_i)^n  in PSD order, for projection matrices.
Exact checks (not a proof) of the rank-one odd-n counterexample family and its parameter transition.
  A. n = 3, 5, 7: v_i = (n-1) e_i - z, P_i = v_i v_i^T / d, d = n^2 - 3n + 3. Each P_i is a symmetric idempotent of
     rank one; the symmetrized product G is computed by enumerating all n! words, exactly over Q; the quadratic form
     z^T (M^n - G) z equals (1 - (n(n-2))^(n-1)) / (n^(n-1) d^n) < 0.
  B. The equiangular family with Gram (1-t)I + tJ: numerically, M^n - G fails to be PSD exactly for
     -1/(n-1) < t < -1/(2n-1) (odd n = 3, 5), and is PSD outside it; even n = 2, 4 show no failure on the grid.
Exit 0 and final line `bad= 0` when every check agrees."""
import itertools, sys
from fractions import Fraction as Fr
import numpy as np

def matmul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

bad = 0
for n in (3, 5, 7):
    d = n * n - 3 * n + 3
    z = [Fr(1)] * n
    v = [[Fr(n - 1) * (1 if j == i else 0) - 1 for j in range(n)] for i in range(n)]
    P = [[[v[i][a] * v[i][b] / d for b in range(n)] for a in range(n)] for i in range(n)]
    for i in range(n):
        if matmul(P[i], P[i]) != P[i] or any(P[i][a][b] != P[i][b][a] for a in range(n) for b in range(n)):
            bad += 1; print('not a projection', n, i)
    M = [[sum(P[i][a][b] for i in range(n)) / n for b in range(n)] for a in range(n)]
    Mn = M
    for _ in range(n - 1):
        Mn = matmul(Mn, M)
    # G via rank-one words: P_s(n)...P_s(1) = c * v_s(n) v_s(1)^T with c = prod of adjacent inner products / d^n
    G = [[Fr(0)] * n for _ in range(n)]
    ip = [[sum(v[i][k] * v[j][k] for k in range(n)) for j in range(n)] for i in range(n)]
    count = 0
    for s in itertools.permutations(range(n)):
        c = Fr(1)
        for a, b in zip(s[1:], s[:-1]):
            c *= ip[a][b]
        c /= Fr(d) ** n
        last, first = s[-1], s[0]
        for a in range(n):
            for b in range(n):
                G[a][b] += c * v[last][a] * v[first][b]
        count += 1
    G = [[x / count for x in row] for row in G]
    if n == 3:  # direct word products as an independent check of the rank-one formula
        G2 = [[Fr(0)] * n for _ in range(n)]
        for s in itertools.permutations(range(n)):
            W = P[s[0]]
            for k in s[1:]:
                W = matmul(P[k], W)
            for a in range(n):
                for b in range(n):
                    G2[a][b] += W[a][b] / 6
        if G2 != G:
            bad += 1; print('rank-one word formula mismatch')
    D = [[Mn[a][b] - G[a][b] for b in range(n)] for a in range(n)]
    q = sum(z[a] * D[a][b] * z[b] for a in range(n) for b in range(n))
    closed = Fr(1 - (n * (n - 2)) ** (n - 1), n ** (n - 1) * d ** n)
    print(f'n={n} d={d} zT(M^n-G)z = {q} closed form {closed} words {count}')
    if q != closed or q >= 0:
        bad += 1; print('A mismatch', n)

def family_min_eig(n, t):
    Gram = (1 - t) * np.eye(n) + t * np.ones((n, n))
    L = np.linalg.cholesky(Gram)          # rows are unit vectors u_i with the prescribed Gram matrix
    U = L
    P = [np.outer(U[i], U[i]) for i in range(n)]
    M = sum(P) / n
    Mn = np.linalg.matrix_power(M, n)
    G = np.zeros((n, n))
    for s in itertools.permutations(range(n)):
        W = np.eye(n)
        for k in s:
            W = P[k] @ W
        G += W
    G /= len(list(itertools.permutations(range(n))))
    return np.linalg.eigvalsh((Mn - G + (Mn - G).T) / 2).min()

for n in (2, 3, 4, 5):
    lo, hi = -1 / (n - 1), -1 / (2 * n - 1)
    grid = [lo + (1 - lo) * k / 400 for k in range(1, 400)]
    fails = [t for t in grid if family_min_eig(n, t) < -1e-12]
    if n % 2 == 1:
        ok = all(lo < t < hi + 1e-9 for t in fails) and any(t < hi - 0.01 for t in fails)
        inside = [t for t in grid if lo + 1e-3 < t < hi - 1e-3]
        ok = ok and all(t in fails for t in inside)
    else:
        ok = not fails
    print(f'n={n} family: failing t in ({min(fails) if fails else None}, {max(fails) if fails else None}); predicted ({lo:.4f}, {hi:.4f}) for odd n; {"ok" if ok else "BAD"}')
    if not ok:
        bad += 1
print('bad=', bad)
sys.exit(0 if bad == 0 else 1)
