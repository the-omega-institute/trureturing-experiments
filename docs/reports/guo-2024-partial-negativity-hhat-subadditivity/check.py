#!/usr/bin/env python3
"""Y. Guo, arXiv:2301.00334v3 (Results Phys. 57 (2024) 107430), conjecture (24): hhat(rho_AB) <= hhat(rho_A) + hhat(rho_B)
for every bipartite state, where hhat(X) = sqrt(delta1 delta2) with delta1 >= delta2 the two largest eigenvalues.
Exact check (Fractions) of the family rho(d, t) = ((1 - t)/2)(P0 + P1) + (t/d^2) I on C^d (x) C^d, P_r the projector onto
Omega_r = d^(-1/2) sum_i |i, i + r mod d>, for d = 5..12 and t in {0, 1/10}:
  (P0 + P1)^2 = P0 + P1 (so P0 + P1 is a rank-2 projector, trace 2), both marginals are I/d,
  hence delta1 = delta2 = (1 - t)/2 + t/d^2 > 1/d = delta1(I/d) = delta2(I/d), and the gap a - 2/d > 0 when t < d(d-4)/(d^2-2).
Pure Python. usage: check.py"""
from fractions import Fraction as Q
import sys

def family(d, t):
    n = d * d; idx = lambda i, j: i * d + j
    P = [[Q(0)] * n for _ in range(n)]
    for r in (0, 1):
        for i in range(d):
            for k in range(d):
                P[idx(i, (i + r) % d)][idx(k, (k + r) % d)] += Q(1, d)
    rho = [[(1 - t) / 2 * P[a][b] + (t / (d * d) if a == b else 0) for b in range(n)] for a in range(n)]
    return P, rho

def matmul(A, B):
    n = len(A); return [[sum(A[i][k] * B[k][j] for k in range(n) if A[i][k]) for j in range(n)] for i in range(n)]

def marginals(rho, d):
    A = [[sum(rho[i * d + j][k * d + j] for j in range(d)) for k in range(d)] for i in range(d)]
    B = [[sum(rho[i * d + j][i * d + l] for i in range(d)) for l in range(d)] for j in range(d)]
    return A, B

ok = True
for d in range(5, 13):
    for t in (Q(0), Q(1, 10)):
        P, rho = family(d, t)
        proj = matmul(P, P) == P
        tr = sum(P[a][a] for a in range(d * d))
        A, B = marginals(rho, d)
        eye = [[Q(1, d) if i == k else Q(0) for k in range(d)] for i in range(d)]
        a = (1 - t) / 2 + t / (d * d)
        gap = a - Q(2, d)
        good = proj and tr == 2 and A == eye and B == eye and gap > 0
        ok &= good
        print(f"d={d} t={t}: (P0+P1)^2=P0+P1 {proj}, tr(P0+P1)={tr}, marginals I/d {A == eye and B == eye}, "
              f"hhat(rho)={a} > hhat(rho_A)+hhat(rho_B)={Q(2, d)}, gap {gap}: {'PASS' if good else 'FAIL'}")
print("RESULT PASS" if ok else "RESULT FAIL"); print("EXIT=0" if ok else "EXIT=1"); sys.exit(0 if ok else 1)
