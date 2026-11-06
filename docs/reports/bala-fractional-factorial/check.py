#!/usr/bin/env python3
"""Exact transcription and finite checks of Bala's three Chebyshev samplings.

Sources: https://oeis.org/A276098/a276098.pdf, page 5;
https://oeis.org/A276100/internal through A276102/internal.
Gamma factors are formal symbols reduced by Gamma(x+1) = x Gamma(x).
No numerical approximation of Gamma, pi, or any other real number is used.
Finite checks do not establish the universally quantified conjectures.
"""

from fractions import Fraction
from hashlib import sha256
from math import factorial, gcd, prod


FAMILIES = (
    ("A276100", 2, 15, 5, 3, 7, (1, 114688, 77636318760,
        62505037015810048, 53837289804317953893960,
        48066503353826060675410034688)),
    ("A276101", 3, 10, 5, 2, 3, (1, 1458, 9723402, 77636318760,
        665145965903562, 5915482311008318958)),
    ("A276102", 5, 6, 3, 2, 1, (1, 50, 8250, 1636250, 349456250,
        77636318760)),
)


def gamma_normal_form(x):
    """Return c,r representing Gamma(1+x) = c Gamma(1+r), 0 <= r < 1."""
    assert isinstance(x, Fraction) and x >= 0
    whole = x.numerator // x.denominator
    residual = x - whole
    den = residual.denominator
    coefficient = Fraction(
        prod(residual.numerator + den * j for j in range(1, whole + 1)),
        den ** whole,
    )
    return coefficient, residual


def gamma_value(n, d, a, b, c):
    """Evaluate the literal Gamma ratio, canceling identical residual symbols."""
    lo, lo_symbol = gamma_normal_form(Fraction(n, d))
    hi, hi_symbol = gamma_normal_form(Fraction(a * n, d))
    assert lo_symbol == hi_symbol, (n, d, lo_symbol, hi_symbol)
    return Fraction(factorial(a * n), factorial(b * n) * factorial(c * n)) * lo / hi


def progression_value(n, d, a, b, c, length):
    return Fraction(
        d ** (length * n) * factorial(a * n),
        factorial(b * n) * factorial(c * n)
        * prod(n + d * k for k in range(1, length * n + 1)),
    )


def ordinary_value(m):
    return Fraction(factorial(30 * m) * factorial(m),
        factorial(15 * m) * factorial(10 * m) * factorial(6 * m))


def delta(x):
    return sum(sign * ((slope * x).numerator // (slope * x).denominator)
        for sign, slope in ((1, 30), (1, 1), (-1, 15), (-1, 10), (-1, 6)))


def main():
    for name, d, a, b, c, length, source_terms in FAMILIES:
        assert a == 30 // d and a == 1 + d * length and a - b - c == length
        terms = []
        for n in range(401):
            value = gamma_value(n, d, a, b, c)
            assert value == progression_value(n, d, a, b, c, length), (name, n)
            assert value > 0 and value.denominator == 1, (name, n, value)
            if n % d == 0:
                assert value == ordinary_value(n // d), (name, n, "ordinary branch")
            terms.append(value.numerator)
        assert tuple(terms[:len(source_terms)]) == source_terms, (name, "source fidelity")
        digest = sha256(("\n".join(map(str, terms)) + "\n").encode()).hexdigest()
        print(f"{name}: source first {len(source_terms)} terms reproduced exactly")
        print(f"{name}: initial values = {', '.join(map(str, terms[:6]))}")
        print(f"{name}: 401 exact Gamma-recurrence/progression comparisons, n=0..400; all positive integers")
        print(f"{name}: values n=0..400 sha256 = {digest}")

    cells = 0
    for j in range(30):
        assert j - j // 2 - j // 3 - j // 5 >= 0
        # Breakpoints and interior representatives of all 30 cells.
        assert delta(Fraction(j, 30)) >= 0
        assert delta(Fraction(2 * j + 1, 60)) >= 0

    for name, d, a, b, c, length, _ in FAMILIES:
        for n in range(81):
            for q in range(2, 101):
                if gcd(q, d) != 1:
                    continue
                count = sum((n + d * k) % q == 0 for k in range(1, length * n + 1))
                z = n * pow(q, -1, d) % d
                x = Fraction(n - q * z, d * q)
                defect = a * n // q - b * n // q - c * n // q - count
                assert defect == delta(x) >= 0, (name, n, q, defect, delta(x))
                cells += 1
    print(f"Shared translated-floor identity: {cells} exact cases; n=0..80, q=2..100 coprime to d")

    # An independently published theorem used as an additional fidelity check:
    # Radcliffe, July 19, 2026, equation (1), OEIS A364175 (a different sequence).
    radcliffe_terms = (1, 36, 3564, 408408, 49697388, 6249195036)
    for n, expected in enumerate(radcliffe_terms):
        lo, r1 = gamma_normal_form(Fraction(2 * n, 3))
        hi, r2 = gamma_normal_form(Fraction(5 * n, 3))
        assert r1 == r2
        assert Fraction(factorial(6 * n), factorial(3 * n) * factorial(2 * n)) * lo / hi == expected
    print("Radcliffe A364175: all six displayed source terms reproduced exactly")
    print("Finite evidence only: universal integrality and literature priority are not certified")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
