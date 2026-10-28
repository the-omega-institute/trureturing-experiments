#!/usr/bin/env python3
"""Exact checks for Englander's A081215 Conjecture XX; no floating point."""

from fractions import Fraction
from math import comb


PREFIX = (
    1, 0, 1, 5, 41, 434, 5713, 90075, 1657009, 34867844,
    826446281, 21794641505, 633095889817, 20088655029078,
    691413758034721, 25657845139503479, 1021273028302258913,
    43404581642184336392, 1961870762757168078553,
)


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def source_a(n):
    """The source's rational expression; check exactness before integer conversion."""
    require(isinstance(n, int) and n >= 0, ("index", n))
    value = Fraction(n ** (n + 1) + (-1) ** n, (n + 1) ** 2)
    require(value.denominator == 1, ("integrality", n))
    return value.numerator


def multiply_coefficients(left, right):
    answer = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            answer[i + j] += x * y
    return answer


def main():
    require(tuple(source_a(n) for n in range(len(PREFIX))) == PREFIX,
            "OEIS A081215 displayed prefix, indices 0..18")
    print("FIDELITY: all 19 displayed OEIS values (n=0..18) agree")

    examples = ((11, 21794641505, 125, 5), (13, 20088655029078, 216, 6))
    for n, stated_value, modulus, remainder in examples:
        value = source_a(n)
        require(value == stated_value, ("source example value", n))
        require(value % modulus == remainder, ("source example remainder", n))
        print(f"FIDELITY: a({n})={value}; a({n}) mod {modulus}={remainder}")

    # Theorem I is a proved source statement, separate from Conjecture XX.
    for n in range(1, 2002):
        require((source_a(n) - (-1) ** n) % n == 0, ("source Theorem I", n))
    print("FIDELITY: source Theorem I holds for all n=1..2001")

    for m in range(2, 1001):
        value = source_a(2 * m + 1)
        require(value % (m ** 3) == m, ("Conjecture XX", m))
        denominator = 4 * (m + 1) ** 2
        # Independent modular-power route, without modular denominator inversion.
        numerator_residue = (pow(2 * m + 1, 2 * m + 2,
                                 denominator * m ** 3) - 1) % (denominator * m ** 3)
        require(numerator_residue == denominator * m, ("modular-power route", m))
    print("SETTLEMENT CHECK: all m=2..1000 pass exact quotient and modular-power routes")
    print("DOMAIN CHECK: 500 even and 499 odd parameters included; no primality restriction")

    for m in range(2, 129):
        h = 2 * m + 2
        value = source_a(2 * m + 1)
        expansion_at_minus_one = sum(comb(h, j) * h ** j * (-1) ** (h - j)
                                     for j in range(h + 1))
        require(expansion_at_minus_one == (h - 1) ** h, ("binomial integrality", m))
        require((expansion_at_minus_one - 1) % (h ** 2) == 0,
                ("denominator divisibility", m))
        tail = sum(comb(h, j) * 2 ** (j - 3) * m ** (j - 3)
                   for j in range(3, h + 1))
        require((1 + 2 * m) ** h - 1 ==
                4 * m * (m + 1) + 4 * m ** 2 * (m + 1) * (2 * m + 1)
                + 8 * m ** 3 * tail, ("cubic truncation", m))
        require((m + 1) ** 2 * (value - m) == 2 * m ** 3 * (m + 1 + tail),
                ("cancel integer factor four", m))
        inverse = 1 - 2 * m + 3 * m ** 2
        require((m + 1) ** 2 * inverse == 1 + m ** 3 * (4 + 3 * m),
                ("inverse polynomial", m))
        require(value - m == m ** 3 * (
            2 * (m + 1 + tail) * inverse - (4 + 3 * m) * (value - m)),
                ("explicit divisibility witness", m))
    print("OUTLINE CHECK: integrality expansion, cubic truncation, integer cancellation,")
    print("and explicit divisibility witness hold for every m=2..128")

    require(multiply_coefficients([1, 2, 1], [1, -2, 3]) == [1, 0, 0, 4, 3],
            "exact polynomial coefficient identity")
    print("SYMBOLIC CHECK: (X+1)^2*(1-2*X+3*X^2)=1+4*X^3+3*X^4")
    print("BOUNDARY: m=1 gives remainder 0, so the source lower bound m>=2 matters")
    require(source_a(3) % 1 == 0, "excluded boundary")
    print("SCOPE: finite checks and one polynomial identity; no Lean proof of the universal claim")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
