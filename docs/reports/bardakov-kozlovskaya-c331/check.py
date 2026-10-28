#!/usr/bin/env python3
"""Exact checks for Bardakov--Kozlovskaya, Conjecture 3.3(1).

The source defines p_n as the product over the n values obtained from the
n-valued multiplication.  The resultant below is the same product, with
auxiliary nth roots a^n = (-1)^n x and b^n = (-1)^n y eliminated.  All
coefficients are integers; no numerical approximation is used.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import sympy as sp

X, Y, Z, A, B, T = sp.symbols("x y z a b t")
E1, E2, E3 = sp.symbols("e1 e2 e3")


def source_p(n: int) -> sp.Poly:
    """Return p_n(z; x, y) from the exact root product/resultant."""
    resultant = sp.resultant(T**n - 1, Z - (A + T * B) ** n, T)
    expanded = sp.Poly(sp.expand(resultant), A, B, Z)
    sign = (-1) ** n
    expression = 0
    for (a_exp, b_exp, z_exp), coefficient in expanded.terms():
        assert a_exp % n == 0 and b_exp % n == 0
        expression += coefficient * (sign * X) ** (a_exp // n) * (
            sign * Y
        ) ** (b_exp // n) * Z**z_exp
    return sp.Poly(sp.expand(expression), X, Y, Z)


def elementary_form(p: sp.Poly) -> sp.Poly:
    """Express the symmetric p in the elementary symmetric variables."""
    expression, remainder, substitutions = sp.symmetrize(
        p.as_expr(), [X, Y, Z], formal=True
    )
    assert remainder == 0
    assert [rhs for _, rhs in substitutions] == [
        X + Y + Z,
        X * Y + X * Z + Y * Z,
        X * Y * Z,
    ]
    names = {str(symbol): symbol for symbol, _ in substitutions}
    return sp.Poly(expression.subs({names["s1"]: E1,
                                     names["s2"]: E2,
                                     names["s3"]: E3}), E1, E2, E3)


def assert_fidelity() -> None:
    expected = {
        4: E1**4 - 2**3 * E1**2 * E2 + 2**4 * E2**2 - 2**7 * E1 * E3,
        6: (
            E1**6
            - 2**2 * 3 * E1**4 * E2
            + 2**4 * 3 * E1**2 * E2**2
            - 2**6 * E2**3
            - 2 * 3**4 * 17 * E1**3 * E3
            - 2**3 * 3**4 * 19 * E1 * E2 * E3
            + 3**3 * 19**3 * E3**2
        ),
        7: (
            E1**7
            - 5 * 7**4 * E1**4 * E3
            + 2 * 7**6 * E1**2 * E2 * E3
            - 7**7 * E2**2 * E3
            + 7**8 * E1 * E3**2
        ),
    }
    for n, expected_expression in expected.items():
        actual = elementary_form(source_p(n))
        assert actual == sp.Poly(sp.expand(expected_expression), E1, E2, E3)
        print(f"fidelity n={n}: source coefficient table reproduced")


def assert_prime_power_congruences() -> None:
    cases = [(2, 1), (3, 1), (2, 2), (5, 1), (2, 3), (7, 1), (3, 2)]
    for prime, exponent in cases:
        n = prime**exponent
        p = source_p(n)
        difference = sp.Poly(
            sp.expand(p.as_expr() - (X + Y + Z) ** n), X, Y, Z
        )
        assert all(int(coefficient) % prime == 0 for coefficient in difference.coeffs())
        elementary = elementary_form(p)
        target = sp.Poly(E1**n, E1, E2, E3)
        delta = elementary - target
        for monomial, coefficient in delta.terms():
            assert int(coefficient) % prime == 0, (n, monomial, coefficient)
        print(f"prime-power n={n}=({prime}^{exponent}): p_n ≡ e1^{n} (mod {prime})")


def assert_non_prime_power_control() -> None:
    """A non-prime-power control prevents the check from hiding the scope."""
    p = source_p(6)
    difference = sp.Poly(sp.expand(p.as_expr() - (X + Y + Z) ** 6), X, Y, Z)
    assert any(int(coefficient) % 2 != 0 for coefficient in difference.coeffs())
    print("control n=6: the unrestricted integer statement is not accepted")


if __name__ == "__main__":
    assert_fidelity()
    assert_prime_power_congruences()
    assert_non_prime_power_control()
    print("ALL_CHECKS_OK")
