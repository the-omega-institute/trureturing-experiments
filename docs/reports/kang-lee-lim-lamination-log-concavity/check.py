#!/usr/bin/env python3
"""Exact checks for Kang--Lee--Lim Conjecture 1 (arXiv:2508.04396v3).

The symbolic mutation check verifies the proposed genus-one extended matrix
specializes to (1+q)^2+q^3.  The coefficient checks use exact integers.
"""
from fractions import Fraction
import sympy as sp

q = sp.Symbol("q")

# Principal matrix and one frozen shear row from the nominated torus witness.
B = sp.Matrix([
    [0, 2, -1, -1],
    [-2, 0, 1, 1],
    [1, -1, 0, 1],
    [1, -1, -1, 0],
])
row = sp.Matrix([[-1, 3, 1, -3]])

# Extended exchange matrix: four mutable rows and one frozen row.
ext = B.col_join(row)


def mutate_matrix(M, k):
    """Fomin--Zelevinsky matrix mutation, exact integer arithmetic."""
    n, m = M.shape
    out = sp.zeros(n, m)
    for i in range(n):
        for j in range(m):
            if i == k or j == k:
                out[i, j] = -M[i, j]
            else:
                out[i, j] = M[i, j] + (abs(M[i, k]) * M[k, j] + M[i, k] * abs(M[k, j])) // 2
    return out

# Two mutations at indices 1 and 2 (zero-based 0 and 1).
ext1 = mutate_matrix(ext, 0)
ext2 = mutate_matrix(ext1, 1)
assert list(ext1[4, :]) == [1, 3, 0, -4]
assert list(ext2[4, :]) == [7, -3, 0, -4]

x1, x2, x3, x4, u = sp.symbols("x1 x2 x3 x4 u")
first = (u * x2**2 + x3 * x4) / x1
# After the first mutation, the column-2 positive/negative products are
# (x1')^2 and u^3*x3*x4, respectively.
second = (first**2 + u**3 * x3 * x4) / x2
specialized = sp.factor(second.subs({x1: 1, x2: 1, x3: 1, x4: 1, u: q}))
assert sp.simplify(specialized - ((1 + q)**2 + q**3)) == 0

# Fidelity: the source's octagon example gives this c-polynomial exactly.
fidelity = sp.Poly(2*q**3 + 4*q**2 + 2*q + 1, q)
assert fidelity.all_coeffs() == [2, 4, 2, 1]
assert [fidelity.nth(i) for i in range(4)] == [1, 2, 4, 2]


def family_coefficients(m):
    p = sp.Poly((1 + q)**2 + q**m, q)
    return [int(p.nth(i)) for i in range(m + 1)]

# The nominated cases m=3 and m=4 fail at their last internal index.
# The draft makes no claim for larger m.
for m in (3, 4):
    coeff = family_coefficients(m)
    i = m - 1
    assert coeff[i] ** 2 < coeff[i - 1] * coeff[i + 1]
    assert all(c >= 0 for c in coeff)

print("source_octagon_coefficients", [fidelity.nth(i) for i in range(4)])
print("mutated_specialization", sp.expand(specialized))
for m in (3, 4):
    coeff = family_coefficients(m)
    i = m - 1
    print(f"m={m} coefficients={coeff} failing_index={i} lhs={coeff[i]**2} rhs={coeff[i-1]*coeff[i+1]}")
print("ALL_CHECKS_OK")
