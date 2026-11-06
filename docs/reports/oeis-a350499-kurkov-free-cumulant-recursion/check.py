#!/usr/bin/env python3
"""Exact symbolic n<=8 check for OEIS A350499 (no floating point)."""
from __future__ import annotations
import sympy as sp

N = 8
MAX = N + 1
t = sp.Symbol("t")
m = {j: sp.Symbol(f"m{j}") for j in range(1, 2 * N + 1)}

Series = dict[int, sp.Expr]
def add(a: Series, b: Series) -> Series:
    return {i: sp.expand(a.get(i, 0) + b.get(i, 0)) for i in range(MAX + 1)}
def mul(a: Series, b: Series) -> Series:
    out = {i: sp.Integer(0) for i in range(MAX + 1)}
    for i, ai in a.items():
        if ai == 0:
            continue
        for j, bj in b.items():
            if i + j <= MAX and bj != 0:
                out[i + j] += ai * bj
    return {i: sp.expand(v) for i, v in out.items()}
def power(a: Series, k: int) -> Series:
    out = {0: sp.Integer(1)}
    for _ in range(k):
        out = mul(out, a)
    return out

def coeff_M_of_U(U: Series, r: int) -> sp.Expr:
    value = U.get(r, 0)
    for j in range(1, r):
        value += m[j] * power(U, j + 1).get(r, 0)
    return sp.expand(value)

# Build M^{-1}(t)=U(t) coefficient by coefficient from M(U)=t.
U: Series = {1: sp.Integer(1)}
for r in range(2, MAX + 1):
    u_r = -coeff_M_of_U(U, r)
    U[r] = sp.expand(u_r)

# C=t/U: invert the unit series U/t coefficientwise.
unit = {q: U.get(q + 1, 0) for q in range(MAX)}
C: Series = {0: sp.Integer(1)}
for q in range(1, N + 1):
    C[q] = sp.expand(-sum(unit[i] * C.get(q - i, 0)
                            for i in range(1, min(q, MAX - 1) + 1)))
source_coeff = {n: C[n] for n in range(1, N + 1)}

# Kurkov's recurrence. Columns through 2N are enough for the first N rows.
R: dict[tuple[int, int], sp.Expr] = {}
for k in range(1, 2 * N + 1):
    R[(1, k)] = m[k]
for n in range(2, N + 1):
    for k in range(1, 2 * N - n + 2):
        R[(n, k)] = sp.expand(
            R[(n - 1, k + 1)]
            - sum(R[(j, k)] * R[(n - j, 1)] for j in range(1, n))
        )

for n in range(1, N + 1):
    assert sp.expand(R[(n, 1)] - source_coeff[n]) == 0, n

# Fidelity against A350499's displayed c_1,...,c_5.
expected = {
    1: m[1],
    2: m[2] - m[1] ** 2,
    3: m[3] - 3 * m[2] * m[1] + 2 * m[1] ** 3,
    4: m[4] - 2 * m[2] ** 2 - 4 * m[3] * m[1] + 10 * m[2] * m[1] ** 2 - 5 * m[1] ** 4,
    5: m[5] - 5 * m[2] * m[3] - 5 * m[4] * m[1] + 15 * m[2] ** 2 * m[1] + 15 * m[3] * m[1] ** 2 - 35 * m[2] * m[1] ** 3 + 14 * m[1] ** 5,
}
for n, value in expected.items():
    assert sp.expand(source_coeff[n] - value) == 0, n

# Independent reversion check through degree 8.
for r in range(1, N + 1):
    target = sp.Integer(1) if r == 1 else sp.Integer(0)
    assert sp.expand(coeff_M_of_U(U, r) - target) == 0, r

print("A350499 exact symbolic recurrence/reversion check")
for n in range(1, N + 1):
    print(f"n={n}: R(n,1) = {sp.expand(R[(n, 1)])}")
print("Verified R(n,1) = [t^n](t/M^(-1)(t)) for every 1 <= n <= 8.")
print("Verified the OEIS examples c_1 through c_5 exactly.")
print("ALL_CHECKS_OK")
