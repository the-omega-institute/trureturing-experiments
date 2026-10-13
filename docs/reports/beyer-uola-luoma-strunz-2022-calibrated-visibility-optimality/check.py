#!/usr/bin/env python3
"""Optimality of the Beyer-Uola-Luoma-Strunz visibility bound under the fluctuation condition.

Source: K. Beyer, R. Uola, K. Luoma, W. T. Strunz, "Joint measurability in nonequilibrium quantum
thermodynamics", arXiv:2111.02854v2, Phys. Rev. E 106, L022101 (2022), p. 4:
"we conjecture that our construction is optimal to satisfy condition (ii) for the given model Eq. (5)".
Model: A_a = lam P_a + (1-lam)/d I, B^U_b = gam U^dag P_b U + (1-gam)/d I; bound Eq. (8)-(9):
gam <= 2 kap / (d + 2 kap - d kap), kap = 2 sqrt(lam + (1-lam)/d) sqrt((1-lam)/d) + (d-2)(1-lam)/d.

Part 1 (sympy, exact): with r = lam + s, s = (1-lam)/d, t = gam/2 + (1-gam)/d,
  kap = 1 - (sqrt r - sqrt s)^2, and  gam/2 <= kap t  <=>  gam (d + 2 kap - d kap) <= 2 kap.
  The left inequality is what every admissible joint POVM forces at the Hadamard-rotated basis
  (diagonal entries fixed by condition (ii) at the basis states, the off-diagonal (0,1) entry summed
  through the second marginal, and the 2x2 principal minor of each PSD effect).
Part 2 (exact rationals): d = 3, lam = 1/2, gam = 4/5 violates the bound (2/5 > 7/18 = kap t).
Part 3 (float64): at gam equal to the bound, for random Haar bases (d = 2..6), the effects
  C_b = (gam/kap) Q_b + u I - (gam (1-kap)/kap) D_b are PSD and W_ab = R_a C_b R_a meets both
  marginals and the diagonal identity (W_ab)_kk = (A_a)_kk (B_b)_kk (the source's own sufficiency).
Prints ALL_OK.
"""
import sys
import sympy as sp
import numpy as np
from fractions import Fraction as Fr

def fail(m):
    print("FAIL", m); sys.exit(1)

# Part 1
lam, gam, d = sp.symbols('lam gam d', positive=True)
s = (1 - lam) / d; r = lam + s
kap = 2 * sp.sqrt(r) * sp.sqrt(s) + (d - 2) * s
if sp.simplify(sp.expand(kap - (1 - (sp.sqrt(r) - sp.sqrt(s))**2))) != 0: fail("kappa identity")
K = sp.Symbol('K', positive=True)
t = gam / 2 + (1 - gam) / d
lhs = sp.expand((K * t - gam / 2) * 2 * d)            # 2d (K t - gam/2) >= 0
rhs = sp.expand(2 * K - gam * (d + 2 * K - d * K))     # 2K - gam (d + 2K - dK) >= 0
if sp.simplify(lhs - rhs) != 0: fail("bound equivalence")
print("part1: kappa = 1 - (sqrt r - sqrt s)^2 and gam/2 <= kap t <=> gam (d + 2 kap - d kap) <= 2 kap (exact)")

# Part 2
d3, l3, g3 = 3, Fr(1, 2), Fr(4, 5)
s3 = (1 - l3) / d3; r3 = l3 + s3
# sqrt(r s) = sqrt(1/9) = 1/3 exactly
assert r3 * s3 == Fr(1, 9)
k3 = 2 * Fr(1, 3) + (d3 - 2) * s3
t3 = g3 / 2 + (1 - g3) / d3
bound = 2 * k3 / (d3 + 2 * k3 - d3 * k3)
if not (k3 == Fr(5, 6) and t3 == Fr(7, 15) and bound == Fr(10, 13)): fail("instance values")
if not (g3 / 2 > k3 * t3 and g3 > bound): fail("instance violation")
print(f"part2: d=3 lam=1/2 gam=4/5: kap={k3}, t={t3}, kap*t={k3*t3} < gam/2={g3/2}; bound={bound}")

# Part 3
rng = np.random.default_rng(20211104)
for dd in range(2, 7):
    for lamv in (0.2, 0.5, 0.8):
        sv = (1 - lamv) / dd; rv = lamv + sv
        kv = 2 * np.sqrt(rv * sv) + (dd - 2) * sv
        gv = 2 * kv / (dd + 2 * kv - dd * kv); uv = (1 - gv) / dd
        A = [np.diag([rv if k == a else sv for k in range(dd)]).astype(complex) for a in range(dd)]
        R = [np.diag([np.sqrt(rv) if k == a else np.sqrt(sv) for k in range(dd)]).astype(complex) for a in range(dd)]
        for _ in range(20):
            Z = rng.normal(size=(dd, dd)) + 1j * rng.normal(size=(dd, dd))
            Q, _r = np.linalg.qr(Z)
            Qs = [np.outer(Q[:, b], Q[:, b].conj()) for b in range(dd)]
            B = [gv * Qs[b] + uv * np.eye(dd) for b in range(dd)]
            C = [(gv / kv) * Qs[b] + uv * np.eye(dd) - (gv * (1 - kv) / kv) * np.diag(np.diag(Qs[b]).real) for b in range(dd)]
            for b in range(dd):
                if np.linalg.eigvalsh(C[b]).min() < -1e-10: fail(f"C not PSD d={dd}")
            W = [[R[a] @ C[b] @ R[a] for b in range(dd)] for a in range(dd)]
            for a in range(dd):
                if np.abs(sum(W[a]) - A[a]).max() > 1e-10: fail("first marginal")
            for b in range(dd):
                if np.abs(sum(W[a][b] for a in range(dd)) - B[b]).max() > 1e-10: fail("second marginal")
                for a in range(dd):
                    if np.abs(np.diag(W[a][b]) - np.diag(A[a]) * np.diag(B[b])).max() > 1e-10: fail("diagonal identity")
    print(f"part3: d={dd}: construction at the bound is a valid calibrated joint POVM for 60 random bases")
print("ALL_OK")
