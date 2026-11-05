#!/usr/bin/env python3
"""Exact counterexample check for the N <= m extension of Proposition 4.16."""
from fractions import Fraction as F
import sympy as sp


def assert_equal(name, lhs, rhs):
    if isinstance(lhs, sp.MatrixBase) or isinstance(rhs, sp.MatrixBase):
        ok = sp.simplify(lhs - rhs) == sp.zeros(*lhs.shape)
    else:
        ok = sp.simplify(lhs - rhs) == 0
    if not ok:
        raise AssertionError(f"{name}: {lhs} != {rhs}")


def assert_true(name, condition):
    if not condition:
        raise AssertionError(name)


# Symbolic family. Rows are spectral modes and columns are sites.
a, b, gamma = sp.symbols("a b gamma", real=True)
A = sp.Matrix([[1, 2], [2, 1]])
pi = sp.Matrix([sp.Rational(1, 2), sp.Rational(1, 2)])
D = sp.Matrix([[a, -a], [b, -b]])
B = A + D
q = sp.Matrix([(2 * b - a) / (a + b), (2 * a - b) / (a + b)])
f = sp.Matrix([1, gamma])

assert_equal("D pi", D * pi, sp.zeros(2, 1))
assert_equal("q sums to one", q[0] + q[1], sp.Integer(1))
assert_equal("Dq proportional to Aq", D * q, (b - a) * A * q)
assert_equal("Bq proportional to Aq", B * q, (1 + b - a) * A * q)

u_pi = (f.T * A * pi)[0]
s_pi = (f.T * (A + B) * pi)[0]
u_q = (f.T * A * q)[0]
s_q = (f.T * (A + B) * q)[0]
assert_equal("C(pi,gamma)", u_pi / s_pi, sp.Rational(1, 2))
assert_equal("C(q,gamma)", u_q / s_q, 1 / (2 + b - a))

# H1: 1 and gamma are linearly independent as functions on (0,1).
assert_true("H1 witness values", sp.Matrix([[1, 1], [sp.Rational(1, 3), sp.Rational(2, 3)]]).det() != 0)
assert_true("A full column rank", A.det() != 0)

# Exact admissible witness a=1/2, b=3/4.
av, bv = F(1, 2), F(3, 4)
Aq = A.subs({a: av, b: bv}) * q.subs({a: av, b: bv})
Bv = B.subs({a: av, b: bv})
qv = q.subs({a: av, b: bv})
Dv = D.subs({a: av, b: bv})
assert_true("a,b are in the strict family", F(0) < av < F(1) and F(0) < bv < F(1) and av / 2 < bv < 2 * av and av != bv)
assert_true("q is an interior distribution", all(x > 0 for x in qv))
assert_true("A and B are entrywise positive", all(x > 0 for x in list(A) + list(Bv)))
assert_true("B1 independence", A.det() != 0)
assert_equal("exact q", qv, sp.Matrix([sp.Rational(4, 5), sp.Rational(1, 5)]))
assert_equal("exact Aq", Aq, sp.Matrix([sp.Rational(6, 5), sp.Rational(9, 5)]))
assert_equal("exact B", Bv, sp.Matrix([[sp.Rational(3, 2), sp.Rational(3, 2)], [sp.Rational(11, 4), sp.Rational(1, 4)]]))
assert_equal("exact Dq", Dv * qv, sp.Rational(1, 4) * Aq)

# B2 fails for every permutation/involution and positive weight: neither B column
# is proportional to either A column.  Determinants test proportionality exactly.
Acols = [A[:, 0], A[:, 1]]
Bcols = [Bv[:, 0], Bv[:, 1]]
for zi, col in enumerate(Bcols):
    for wi, acol in enumerate(Acols):
        assert_true(f"nonduality B column {zi} vs A column {wi}", sp.det(sp.Matrix.hstack(col, acol)) != 0)

# H2 holds on (0,1): all entries of A and B are positive and f=(1,gamma)>0,
# hence every s_z(gamma) is positive.  Check it symbolically at the exact witness.
for z in range(2):
    s_z = sp.expand((f.T * (A + Bv)[:, z])[0])
    assert_true(f"H2 s_{z}", all(coef > 0 for coef in sp.Poly(s_z, gamma).all_coeffs()))

# Two distinct insensitive distributions at levels 1/2 and 4/9 refute uniqueness.
C_pi_exact = (f.T * A * pi)[0] / (f.T * (A + Bv) * pi)[0]
C_q_exact = (f.T * A * qv)[0] / (f.T * (A + Bv) * qv)[0]
assert_equal("exact C(pi,gamma)", C_pi_exact, sp.Rational(1, 2))
assert_equal("exact C(q,gamma)", C_q_exact, sp.Rational(4, 9))
assert_true("distinct levels", sp.Rational(1, 2) != sp.Rational(4, 9))
assert_true("distinct distributions", pi != qv)

# A second exact point confirms the displayed strict inequalities describe a
# nonempty open family in the two-dimensional space D*pi=0.
for av2, bv2 in ((F(1, 3), F(1, 2)), (F(3, 5), F(7, 10))):
    assert_true("second family point admissible", F(0) < av2 < F(1) and F(0) < bv2 < F(1) and av2 / 2 < bv2 < 2 * av2 and av2 != bv2)

print("symbolic_family_identities=ok")
print("exact_witness=a=1/2,b=3/4")
print("C_pi=1/2,C_q=4/9")
print("open_family_constraints=0<a<1,0<b<1,a/2<b<2a,a!=b")
print("ALL_CHECKS_OK")
