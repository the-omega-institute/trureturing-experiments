#!/usr/bin/env python3
"""Exact checks for the W-state Bell-state-measurement candidate.

All matrix entries and inequalities below are evaluated over SymPy's exact
integers, rationals, and radicals; no floating-point values are used.
"""

from fractions import Fraction
from itertools import product

import sympy as sp


I = sp.I
zero = sp.Integer(0)
one = sp.Integer(1)
sqrt = sp.sqrt

pauli = {
    0: sp.eye(2),
    1: sp.Matrix([[0, 1], [1, 0]]),
    2: sp.Matrix([[0, -I], [I, 0]]),
    3: sp.Matrix([[1, 0], [0, -1]]),
}


def kron3(a, b, c):
    return sp.kronecker_product(a, b, c)


def ket_w():
    return sp.Matrix([0, 1, 1, 0, 1, 0, 0, 0]) / sqrt(3)


psi = ket_w()
rho = psi * psi.conjugate().T


def tau(mu, nu, state=rho):
    return sp.Matrix([
        sp.simplify(sp.trace(state * kron3(pauli[mu], pauli[k], pauli[nu])))
        for k in (1, 2, 3)
    ])


def norm_sq(v):
    return sp.simplify((v.conjugate().T * v)[0])


def lambda_value(u, v):
    s = sp.simplify(norm_sq(u) + norm_sq(v))
    g = sp.simplify((u.conjugate().T * v)[0])
    return sp.simplify((s + sqrt(s * s - 4 * g * g)) / 2)


def assert_exact(condition, message):
    if not bool(condition):
        raise AssertionError(message)


# Fidelity with the paper's displayed W correlations.
assert_exact(tau(3, 3) == sp.Matrix([0, 0, -1]), "tau_33 != -z")
assert_exact(tau(1, 3) == sp.Matrix([sp.Rational(2, 3), 0, 0]), "tau_13 != 2x/3")
u = tau(3, 3)
v = tau(1, 3)
assert_exact(lambda_value(u, v) == sp.Rational(13, 9), "attaining Lambda is not 13/9")

# Exact finite helper check: all index pairs in the standard Cartesian frames
# obey the claimed 13/9 bound.  The analytic preregistered proof must still
# cover arbitrary orthonormal frames; this is only a regression check.
vectors = {(mu, nu): tau(mu, nu) for mu in range(4) for nu in range(4)}
for (mu, nu), (mu2, nu2) in product(vectors, repeat=2):
    value = lambda_value(vectors[(mu, nu)], vectors[(mu2, nu2)])
    assert_exact(value <= sp.Rational(13, 9), f"Cartesian pair exceeds 13/9: {(mu, nu), (mu2, nu2)}")

# Exact check of the state-specific quadratic identity on a small rational
# unit-vector family.  T(a,c) is built directly from the W state.
def T(a, c):
    result = sp.zeros(3, 1)
    for k in range(3):
        observable = sum((a[i] * pauli[i + 1] for i in range(3)), sp.zeros(2))
        observable_c = sum((c[i] * pauli[i + 1] for i in range(3)), sp.zeros(2))
        result[k] = sp.simplify(sp.trace(rho * kron3(observable, pauli[k + 1], observable_c)))
    return result


unit_vectors = [
    (sp.Integer(1), 0, 0),
    (0, sp.Integer(1), 0),
    (0, 0, sp.Integer(1)),
    (sp.Rational(3, 5), sp.Rational(4, 5), 0),
    (sp.Rational(3, 5), 0, sp.Rational(4, 5)),
    (0, sp.Rational(3, 5), sp.Rational(4, 5)),
]
D = sp.diag(2, 2, -1)
for a, c in product(unit_vectors, repeat=2):
    av = sp.Matrix(a)
    cv = sp.Matrix(c)
    lhs = sp.simplify(9 * norm_sq(T(av, cv)))
    rhs = sp.simplify((av.T * D * cv)[0] ** 2 + 4 * av[2] ** 2 + 4 * cv[2] ** 2)
    assert_exact(lhs == rhs, f"quadratic identity failed for {a}, {c}")

# Traceless middle Pauli factors make white-noise contributions vanish, so the
# displayed correlations scale exactly by v.
v_symbol = sp.symbols("v", real=True)
rho_v = v_symbol * rho + (1 - v_symbol) * sp.eye(8) / 8
for mu, nu in ((3, 3), (1, 3), (0, 0), (2, 1)):
    assert_exact(tau(mu, nu, rho_v) == v_symbol * tau(mu, nu), f"visibility scaling failed for {(mu, nu)}")

print("W correlations: tau_33 = (0, 0, -1), tau_13 = (2/3, 0, 0)")
print("Cartesian index-pair bound and quadratic identity checks: OK")
print("White-noise visibility scaling: OK")
print("ALL_CHECKS_OK")
