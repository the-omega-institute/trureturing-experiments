#!/usr/bin/env python3
"""Exact fidelity and common-metric checks for arXiv:2606.12865v1, Q6.

Requires SymPy. All evidence is symbolic or rational; no floating point is used.
This is an experiment, not a Lean proof or a prior-art determination.
"""

import platform
from fractions import Fraction

import sympy as s


def eq(actual, expected):
    """Compare every entry by exact simplification."""
    difference = actual - expected
    assert all(s.simplify(s.expand(x)) == 0 for x in difference), difference


def beta(T, m):
    return sum(
        ((-1) ** (m - k) * s.binomial(m, k) * T.H**k * T**k
         for k in range(m + 1)), s.zeros(T.rows)
    )


def canonical_inverse(T, m):
    return sum(
        ((-1) ** (k - 1) * s.binomial(m, k) * T.H**k * T**(k - 1)
         for k in range(1, m + 1)), s.zeros(T.rows)
    )


def candidate(T):
    return 2 * T.H - T.H**2 * T


def reflexive(T, S):
    P, Q = T * S, S * T
    eq(T * S * T, T)
    eq(S * T * S, S)
    eq(P**2, P)
    eq(Q**2, Q)
    return P.applyfunc(s.simplify), Q.applyfunc(s.simplify)


print("Python", platform.python_version(), "SymPy", s.__version__)
I, Z = s.eye(2), s.zeros(2)

# Fidelity: apply the source's numbered Lemma 2.8 to its own Q6 example.
phi = (1 + s.sqrt(5)) / 2
a0 = s.sqrt(phi)
T0 = s.Matrix([[a0, 0], [1, 0]])
S0 = canonical_inverse(T0, 2)
eq(S0, candidate(T0))
eq(beta(T0, 2), s.diag(0, 1))
eq(T0 * beta(T0, 2), Z)
P0, Q0 = reflexive(T0, S0)
eq(Q0, I - beta(T0, 2))
eq(Q0, s.diag(1, 0))
eq(P0, s.Matrix([[-1, 2 * a0], [a0 * (1 - s.sqrt(5)) / 2, 2]]))
assert s.simplify(P0[0, 1] - s.conjugate(P0[1, 0])) != 0
print("FIDELITY Lemma 2.8: canonical S, TST=T, STS=S, ST=I-beta_2(T)")
print("FIDELITY source Q6 example: displayed TS reproduced exactly")

# The rational witness is additionally checked using fractions.Fraction only.
def mm(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(2)), Fraction(0))
             for j in range(2)] for i in range(2)]


def transpose(A):
    return [list(row) for row in zip(*A)]


F = Fraction
Tf = [[F(7, 5), F(0)], [F(24, 5), F(0)]]
Af = transpose(Tf)
Af2Tf = mm(mm(Af, Af), Tf)
Sf = [[2 * Af[i][j] - Af2Tf[i][j] for j in range(2)] for i in range(2)]
Pf, Qf = mm(Tf, Sf), mm(Sf, Tf)
assert Sf == [[F(-161, 5), F(48, 5)], [F(0), F(0)]]
assert Pf == [[F(-1127, 25), F(336, 25)], [F(-3864, 25), F(1152, 25)]]
assert Qf == [[F(1), F(0)], [F(0), F(0)]]
assert mm(Pf, Tf) == Tf and mm(Qf, Sf) == Sf
assert mm(Pf, Pf) == Pf and mm(Qf, Qf) == Qf
T = s.Matrix(Tf)
S = candidate(T)
eq(T.H * T, s.diag(25, 0))
eq(T**2, s.Matrix([[s.Rational(49, 25), 0], [s.Rational(168, 25), 0]]))
eq(T.H**2 * T**2, s.diag(49, 0))
eq(beta(T, 2), s.diag(0, 1))
eq(T * beta(T, 2), Z)
eq(T * beta(T, 1), 24 * T)
eq(S, s.Matrix(Sf))
P, Q = reflexive(T, S)
assert T.rank() == 1 and T != Z
assert T * T.H != T.H * T
print("RATIONAL T =", T.tolist())
print("RATIONAL S =", S.tolist())
print("RATIONAL TS =", P.tolist())
print("RATIONAL ST =", Q.tolist())
print("RATIONAL beta_2=diag(0,1); T beta_1=24T; rank=1; nonnormal")

# All real r >= 2, represented exactly by r=2+t with t >= 0.
t = s.symbols("t", nonnegative=True)
r = 2 + t
a = s.sqrt(2 * r - 1) / s.sqrt(r)
assert s.simplify(a - s.sqrt(2 - 1 / r)) == 0
b = (r - 1) / s.sqrt(r)
Tr = s.Matrix([[a, 0], [b, 0]])
Sr = candidate(Tr)
eq(s.Matrix([a**2 + b**2]), s.Matrix([r]))
assert a.is_positive is True and b.is_positive is True
eq(beta(Tr, 2), s.diag(0, 1))
eq(Tr * beta(Tr, 2), Z)
eq(Tr * beta(Tr, 1), (r - 1) * Tr)
eq(Sr, s.Matrix([[a * (2 - r), 2 * b], [0, 0]]))
Pr, Qr = reflexive(Tr, Sr)
eq(Qr, s.diag(1, 0))
eq(s.Matrix([Pr[0, 1], Pr[1, 0]]), s.Matrix([2 * a * b, a * b * (2 - r)]))
eq(Tr.subs(t, 23), T)
assert s.simplify(a * b).is_positive is True
print("FAMILY r=2+t, t>=0: a>0, b>0, a^2+b^2=r")
print("FAMILY beta_2=diag(0,1); T beta_2=0; T beta_1=(r-1)T != 0")
print("FAMILY TST=T; STS=S; P^2=P; Q^2=Q; Q=diag(1,0)")
print("FAMILY P_12=2ab>0; P_21=ab(2-r)<=0; r=25 is the rational witness")

# Every complex inner product in dimension two has such a Hermitian Gram matrix.
p, q = s.symbols("p q", positive=True)
u, v = s.symbols("u v", real=True)
G = s.Matrix([[p, u + s.I * v], [u - s.I * v, q]])
q_residual = G * Qr - Qr.H * G
eq(q_residual, s.Matrix([[0, -u - s.I * v], [u - s.I * v, 0]]))
assert s.solve([s.re(q_residual[0, 1]), s.im(q_residual[0, 1])], [u, v]) == {u: 0, v: 0}
D = G.subs({u: 0, v: 0})
obstruction = s.simplify((D * Pr - Pr.H * D)[0, 1])
eq(s.Matrix([obstruction]), s.Matrix([a * b * (2 * p + t * q)]))
assert (2 * p + t * q).is_positive is True
assert obstruction.is_positive is True
rational_obstruction = (D * P - P.H * D)[0, 1]
eq(s.Matrix([rational_obstruction]), s.Matrix([(336 * p + 3864 * q) / 25]))
print("COMMON METRIC GQ=Q*G forces u=v=0 and positivity forces p,q>0")
print("COMMON METRIC (GP-P*G)_12=ab(2p+tq)>0 for every allowed metric")
print("RATIONAL COMMON METRIC 336p=-3864q is impossible for p,q>0")

# Each idempotent can be orthogonalized individually, including the boundary t=0.
# For any real idempotent R, W=R*R+(I-R)*(I-R) is positive definite:
# x*Wx=||Rx||^2+||(I-R)x||^2, and the two zero vectors imply x=0.
W = (Pr.H * Pr + (I - Pr).H * (I - Pr)).applyfunc(s.simplify)
eq(W * Pr, Pr.H * W)
eq(Qr.H, Qr)
print("INDIVIDUAL METRICS Q is standard-selfadjoint; W=P*P+(I-P)*(I-P) symmetrizes P")
print("INDIVIDUAL METRICS x*Wx=||Px||^2+||(I-P)x||^2>0 for x!=0")
print("ALL_CHECKS_OK")
