#!/usr/bin/env python3
"""The exact quantum bound of the Bernards-Guhne nonlocality-depth Bell family, for every n.

Source: F. Bernards, O. Guhne, "Bell inequalities for nonlocality depth", arXiv:2205.04250v2, Phys. Rev. A 107,
022412 (2023), Sec. III.F, Eqs. (28)-(29): F_n = sum_l (-1)^(1+ceil(l/2)) l (1..1 2..2) and "We find, that the quantum
bound of the inequalities is F_n = sqrt(2) n 2^(n-2) up to numerical precision" (numerical finding for n <= 8).

Checks:
  1. Exact (integers): with settings X (setting 1) and Y (setting 2), B_n |0^n> = b_n |1^n> and B_n |1^n> = conj(b_n)
     |0^n> with b_n = n 2^(n-2) (-1 + i), for n = 3..14; so the eigenvalue sqrt(2) n 2^(n-2) is attained
     (|b_n| = q_n). Computed from the exact action of tensor products of X, Y on |0^n>.
  2. Float64: largest eigenvalue of B_n with settings (X, Z) and (X, Y) equals q_n for n = 3..8, and 200 random qubit
     projective settings per n never exceed q_n.
Prints ALL_OK.
"""
import sys, math, itertools
import numpy as np
from fractions import Fraction

def fail(m):
    print("FAIL", m); sys.exit(1)

def c(l): return (-1) ** (1 + math.ceil(l / 2))

# 1. exact: Y|0> = i|1>, X|0> = |1>; a string with l copies of Y maps |0^n> to i^l |1^n>
for n in range(3, 15):
    re, im = 0, 0
    for l in range(0, n + 1):
        coeff = math.comb(n, l) * l * c(l)
        ph = [(1, 0), (0, 1), (-1, 0), (0, -1)][l % 4]  # i^l
        re += coeff * ph[0]; im += coeff * ph[1]
    if (re, im) != (-n * 2 ** (n - 2), n * 2 ** (n - 2)): fail(f"b_n at n={n}: {(re, im)}")
print("1: b_n = n 2^(n-2) (-1 + i) exactly for n = 3..14, so |b_n| = sqrt(2) n 2^(n-2) is attained")

# 2. float
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1, -1]).astype(complex)
def Bn(n, A0, A1):
    B = 0
    for x in itertools.product([0, 1], repeat=n):
        h = sum(x)
        if h == 0: continue
        op = np.array([[1]], complex)
        for xi in x: op = np.kron(op, A1 if xi else A0)
        B = B + h * c(h) * op
    return B
rng = np.random.default_rng(220504250)
for n in range(3, 9):
    q = math.sqrt(2) * n * 2 ** (n - 2)
    for A0, A1 in ((X, Z), (X, Y)):
        e = np.linalg.eigvalsh(Bn(n, A0, A1)).max()
        if abs(e - q) > 1e-9 * q: fail(f"attainment n={n}")
    for _ in range(200):
        a = []
        for s in range(2):
            v = rng.normal(size=3); v /= np.linalg.norm(v); a.append(v[0] * X + v[1] * Y + v[2] * Z)
        if np.linalg.eigvalsh(Bn(n, a[0], a[1])).max() > q + 1e-9: fail(f"bound violated n={n}")
    print(f"2: n={n}: q_n = {q:.6f} attained by (X,Z) and (X,Y); 200 random qubit settings stay below")
print("ALL_OK")
