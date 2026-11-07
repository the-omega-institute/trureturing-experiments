"""Exact checks for the elliptic-curve candidate from arXiv:2608.24849v1.

The pair representation p(x) + q(x)y, with y^2 = x^3 + ax + b,
keeps all arithmetic in SymPy's exact polynomial ring.  The script checks
the displayed bracket identity, the generating recurrence through a finite
range, and the absence of a single bracket in a bounded exact ansatz.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import sympy as sp


x = sp.symbols("x")


@dataclass(frozen=True)
class Pair:
    p: sp.Expr
    q: sp.Expr


def add(f: Pair, g: Pair) -> Pair:
    return Pair(sp.expand(f.p + g.p), sp.expand(f.q + g.q))


def sub(f: Pair, g: Pair) -> Pair:
    return Pair(sp.expand(f.p - g.p), sp.expand(f.q - g.q))


def mul(f: Pair, g: Pair, P: sp.Expr) -> Pair:
    return Pair(sp.expand(f.p * g.p + f.q * g.q * P),
                sp.expand(f.p * g.q + f.q * g.p))


def derivation(f: Pair, a: int, b: int) -> Pair:
    P = x**3 + a * x + b
    Pp = sp.diff(P, x)
    # D = 2 y d/dx + P'(x) d/dy, reduced modulo y^2=P.
    return Pair(sp.expand(2 * sp.diff(f.q, x) * P + f.q * Pp),
                sp.expand(2 * sp.diff(f.p, x)))


def bracket(f: Pair, g: Pair, a: int, b: int) -> Pair:
    P = x**3 + a * x + b
    return sub(mul(f, derivation(g, a, b), P),
               mul(g, derivation(f, a, b), P))


def assert_pair(actual: Pair, expected: Pair, label: str) -> None:
    if sp.expand(actual.p - expected.p) != 0 or sp.expand(actual.q - expected.q) != 0:
        raise AssertionError(f"{label}: {actual} != {expected}")


def check_curve(a: int, b: int) -> None:
    discriminant = 4 * a**3 + 27 * b**2
    if discriminant == 0:
        raise AssertionError("sample curve is singular")

    one = Pair(sp.Integer(1), sp.Integer(0))
    X = Pair(x, sp.Integer(0))
    Y = Pair(sp.Integer(0), sp.Integer(1))
    P = x**3 + a * x + b

    # Fidelity: the source's formula [yD, x^nD].
    for n in range(2, 8):
        lhs = bracket(Y, Pair(x**n, 0), a, b).p
        rhs = (2 * n - 3) * x**(n + 2) + (2 * n - 1) * a * x**n + 2 * n * b * x**(n - 1)
        if sp.expand(lhs - rhs) != 0:
            raise AssertionError(f"bracket formula failed for n={n}, a={a}, b={b}")

    # The first brackets used by the two-generator induction.
    assert_pair(bracket(one, X, a, b), Pair(2 * Y.p, 2 * Y.q), "[D,xD]")
    assert_pair(bracket(one, Y, a, b), Pair(3 * x**2 + a, 0), "[D,yD]")
    assert_pair(bracket(Y, X, a, b), Pair(-x**3 + a * x + 2 * b, 0), "[yD,xD]")

    # Exact recurrence: x^2, x^3 are obtained from the displayed brackets,
    # and n=2,...,7 solve successively for x^(n+2).
    if sp.expand((bracket(one, Y, a, b).p - a) / 3 - x**2) != 0:
        raise AssertionError("x^2 seed failed")
    if sp.expand(-bracket(Y, X, a, b).p + a * x + 2 * b - x**3) != 0:
        raise AssertionError("x^3 seed failed")
    for n in range(2, 8):
        identity = bracket(Y, Pair(x**n, 0), a, b).p
        solved = sp.expand((identity - (2 * n - 1) * a * x**n - 2 * n * b * x**(n - 1)) / (2 * n - 3))
        if sp.expand(solved - x**(n + 2)) != 0:
            raise AssertionError(f"recurrence failed for n={n}")

    # Remaining basis terms x^k y arise from [D, x^(k+1)D].
    for k in range(0, 8):
        y_coeff = bracket(one, Pair(x ** (k + 1), 0), a, b).q
        if sp.expand(y_coeff - 2 * (k + 1) * x**k) != 0:
            raise AssertionError(f"y-basis induction failed for k={k}")

    # A bounded exact obstruction check.  Unknown coefficients are rationally
    # constrained by a Groebner basis, so no floating-point evidence enters.
    N = 2
    u = sp.symbols(f"u0:{N + 1}")
    v = sp.symbols(f"v0:{N + 1}")
    U = sp.symbols(f"U0:{N + 1}")
    V = sp.symbols(f"V0:{N + 1}")
    f = Pair(sum(u[i] * x**i for i in range(N + 1)),
             sum(v[i] * x**i for i in range(N + 1)))
    g = Pair(sum(U[i] * x**i for i in range(N + 1)),
             sum(V[i] * x**i for i in range(N + 1)))
    B = bracket(f, g, a, b)
    equations = []
    for poly, target in ((B.p, sp.Integer(1)), (B.q, sp.Integer(0))):
        equations.extend(sp.Poly(sp.expand(poly - target), x).all_coeffs())
    variables = list(u) + list(v) + list(U) + list(V)
    basis = sp.groebner(equations, *variables, order="grevlex")
    if any(poly.as_expr() == 1 for poly in basis.polys):
        # No degree-two coefficient ansatz can give the obstruction D.
        return
    raise AssertionError(f"bounded single-bracket obstruction not detected for ({a},{b})")


def main() -> None:
    for a, b in ((0, 1), (1, 1)):
        check_curve(a, b)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
