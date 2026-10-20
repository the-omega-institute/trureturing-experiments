#!/usr/bin/env python3
"""Albright-Hadaway-Pomerance-Jeffries-Lorenzen-Nix, arXiv:2211.01495v2, Conjecture 3.5: the edge derivative of Kemeny's
constant K = sum_{lambda_L != 0} 1/lambda_L (normalized Laplacian of A(t) = A + tB) is smallest on the path P_n with respect
to the non-edge joining its ends. Barbell B_k: two K_k joined end to end by a path through k-1 new vertices (n = 3k - 1),
differentiated at the missing chord between the two attachment vertices.
Part A (floating point, the source's definition): central differences of the eigenvalue sum for k = 3..15.
Part B (exact, Fractions): D_ab K = [sum_j d_j (R_aj + R_bj)]/(2m) - K/m - sum_i d_i (v_i - vbar)^2, with L v = e_a - e_b,
R the effective resistance and K = (1/2m) sum_{i<j} d_i d_j R_ij, for B_8 and B_10 and P_23, P_29.
usage: check.py"""
import sys
from fractions import Fraction as Q
import numpy as np

def barbell(k):
    n = 3 * k - 1; E = set()
    for base in (0, k):
        for i in range(k):
            for j in range(i + 1, k): E.add((base + i, base + j))
    p = [0] + list(range(2 * k, 3 * k - 1)) + [k]
    for i in range(len(p) - 1): E.add(tuple(sorted((p[i], p[i + 1]))))
    return n, sorted(E), 0, k

def path(n): return n, [(i, i + 1) for i in range(n - 1)], 0, n - 1

def kemeny_eig(n, E, a, b, t):
    A = np.zeros((n, n))
    for i, j in E: A[i, j] = A[j, i] = 1.0
    A[a, b] += t; A[b, a] += t
    d = A.sum(1); Dm = np.diag(d ** -0.5); L = np.eye(n) - Dm @ A @ Dm
    w = np.sort(np.linalg.eigvalsh(L)); return float(sum(1 / x for x in w[1:]))

def deriv_eig(n, E, a, b, h=1e-5):
    return (kemeny_eig(n, E, a, b, h) - kemeny_eig(n, E, a, b, -h)) / (2 * h)

def solve(M, rhs):  # exact Gaussian elimination
    n = len(M); A = [row[:] + [r] for row, r in zip(M, rhs)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c] / A[c][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [A[i][n] / A[i][i] for i in range(n)]

def exact_deriv(n, E, a, b):
    deg = [0] * n
    for i, j in E: deg[i] += 1; deg[j] += 1
    m = len(E)
    L = [[Q(0)] * n for _ in range(n)]
    for i, j in E: L[i][i] += 1; L[j][j] += 1; L[i][j] -= 1; L[j][i] -= 1
    # potentials of a unit current from src to snk, with node n-1 grounded
    def potentials(src, snk):
        idx = [i for i in range(n) if i != n - 1]
        M = [[L[i][j] for j in idx] for i in idx]
        rhs = [Q((1 if i == src else 0) - (1 if i == snk else 0)) for i in idx]
        x = solve(M, rhs); v = [Q(0)] * n
        for k, i in enumerate(idx): v[i] = x[k]
        return v
    # effective resistance R_ij = potential difference for a unit current from i to j
    Rm = [[Q(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            v = potentials(i, j); Rm[i][j] = Rm[j][i] = v[i] - v[j]
    K = sum(Q(deg[i] * deg[j]) * Rm[i][j] for i in range(n) for j in range(i + 1, n)) / (2 * m)
    v = potentials(a, b); vbar = sum(deg[i] * v[i] for i in range(n)) / (2 * m)
    var = sum(deg[i] * (v[i] - vbar) ** 2 for i in range(n))
    return sum(deg[j] * (Rm[a][j] + Rm[b][j]) for j in range(n)) / (2 * m) - K / m - var

ok = True
for k in range(3, 16):
    n, E, a, b = barbell(k); nP, EP, aP, bP = path(n)
    dB = deriv_eig(n, E, a, b); dP = deriv_eig(nP, EP, aP, bP)
    closed = -k**4 / 2 + k**3 / 3 + k / 6 + 7 / (6 * k) - 3 / k**2 + 17 / (6 * k**3) - 1 / k**4
    closedP = -((n - 1)**2 - 1)**2 / (6 * (n - 1))
    agree = abs(dB - closed) < 1e-3 * abs(closed) and abs(dP - closedP) < 1e-3 * abs(closedP)
    ok &= agree
    print(f"A k={k} n={n}: D K(B_k) = {dB:.4f} (closed form {closed:.4f}); D_ends K(P_n) = {dP:.4f} (closed form {closedP:.4f}); "
          f"{'counterexample' if dB < dP else 'no violation'}", flush=True)
for k in (8, 10):
    n, E, a, b = barbell(k); nP, EP, aP, bP = path(n)
    dB = exact_deriv(n, E, a, b); dP = exact_deriv(nP, EP, aP, bP)
    closed = Q(-k**4, 2) + Q(k**3, 3) + Q(k, 6) + Q(7, 6 * k) - Q(3, k**2) + Q(17, 6 * k**3) - Q(1, k**4)
    closedP = Q(-((n - 1)**2 - 1)**2, 6 * (n - 1))
    good = dB == closed and dP == closedP and dB < dP; ok &= good
    print(f"B k={k} n={n}: exact D K(B_k) = {dB} = closed form {dB == closed}; exact D_ends K(P_n) = {dP} = closed form "
          f"{dP == closedP}; D K(B_k) < D_ends K(P_n): {dB < dP}", flush=True)
print("RESULT PASS" if ok else "RESULT FAIL"); print("EXIT=0" if ok else "EXIT=1"); sys.exit(0 if ok else 1)
