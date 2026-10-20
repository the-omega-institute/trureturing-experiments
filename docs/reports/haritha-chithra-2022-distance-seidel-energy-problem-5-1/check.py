#!/usr/bin/env python3
"""Haritha-Chithra, arXiv:2210.05940v2, Problem 5.1: E_{D^S}(K_{a,b} - e) > E_{D^S}(K_{a,b}) for a, b >= 2, where
D^S = J - I - 2D and E is the sum of absolute eigenvalues.
Exact part (integers): for 2 <= a <= b <= 12, the deleted-edge matrix equals M - 4(E_uv + E_vu) (only the uv distance
changes, 1 -> 3), and M acts as 3I on the vectors summing to zero on each side (M x = 3x checked on a spanning set).
Floating part: E(K_{a,b}) = 6(n - 2) and E(K_{a,b} - e) > E(K_{a,b}) for 2 <= a <= b <= 12, also under the alternative
reading J - 2D (diagonal 1) of the source's casewise formula. usage: check.py"""
import sys
import numpy as np

def dist(n, E):
    INF = 10**6; D = [[0 if i == j else INF for j in range(n)] for i in range(n)]
    for i, j in E: D[i][j] = D[j][i] = 1
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if D[i][k] + D[k][j] < D[i][j]: D[i][j] = D[i][k] + D[k][j]
    return D

def seidel(n, E, diag=0):
    D = dist(n, E); return [[(1 - 1 if i == j else 1) + (diag if i == j else 0) - 2 * D[i][j] for j in range(n)] for i in range(n)]

ok = True
for a in range(2, 13):
    for b in range(a, 13):
        n = a + b; E = [(i, a + j) for i in range(a) for j in range(b)]
        M = seidel(n, E); N = seidel(n, E[1:]); u, v = E[0]
        diff = [[N[i][j] - M[i][j] for j in range(n)] for i in range(n)]
        want = [[(-4 if {i, j} == {u, v} else 0) for j in range(n)] for i in range(n)]
        exact = diff == want
        # M x = 3x for x = e_i - e_j within one side
        three = all(all(sum(M[r][c] * x[c] for c in range(n)) == 3 * x[r] for r in range(n))
                    for x in ([1 if c == 0 else (-1 if c == i else 0) for c in range(n)] for i in range(1, a)))
        e0 = float(np.abs(np.linalg.eigvalsh(np.array(M, float))).sum()); e1 = float(np.abs(np.linalg.eigvalsh(np.array(N, float))).sum())
        Ma = np.array(seidel(n, E, 1), float); Na = np.array(seidel(n, E[1:], 1), float)
        alt = float(np.abs(np.linalg.eigvalsh(Na)).sum()) > float(np.abs(np.linalg.eigvalsh(Ma)).sum()) + 1e-9
        good = exact and three and abs(e0 - 6 * (n - 2)) < 1e-8 and e1 > e0 + 1e-9 and alt
        ok &= good
        if b == a or a == 2 or not good:
            print(f"a={a} b={b}: N - M = -4(E_uv + E_vu) {exact}; M = 3 on side-difference vectors {three}; E(K) = {e0:.6f} "
                  f"(6(n-2) = {6 * (n - 2)}); E(K - e) = {e1:.6f}; gap {e1 - e0:.6f}; alt reading increases {alt}: {'PASS' if good else 'FAIL'}")
print("RESULT PASS" if ok else "RESULT FAIL"); print("EXIT=0" if ok else "EXIT=1"); sys.exit(0 if ok else 1)
