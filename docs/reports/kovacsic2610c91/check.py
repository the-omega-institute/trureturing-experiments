"""Exact checks for arXiv:2610.14417v1, Proposition 6.2 and Conjecture 9.1."""

from itertools import product

from sympy import Poly, QQ, Rational, cancel, symbols, sympify


t, *variables = symbols("t u1:5")


def rational_parts(value, generators):
    numerator, denominator = cancel(value).as_numer_denom()
    return Poly(numerator, *generators, domain=QQ), Poly(denominator, *generators, domain=QQ)


def valuation(value, generators, weights):
    numerator, denominator = rational_parts(value, generators)
    if numerator.is_zero:
        return None  # The additive valuation is +∞ at zero.

    def polynomial_value(poly):
        return min(
            sum(weight * exponent for weight, exponent in zip(weights, monomial, strict=True))
            for monomial in poly.monoms()
        )

    return polynomial_value(numerator) - polynomial_value(denominator)


def check_kochen(sample, generators, weights):
    k = valuation(sample, generators, weights)
    denominator = cancel(1 + t * sample**2)
    assert denominator != 0
    image = valuation(sample / denominator, generators, weights)
    if k is None:
        assert image is None
        assert valuation(denominator, generators, weights) == 0
        return
    assert valuation(denominator, generators, weights) == min(0, 1 + 2 * k)
    assert image == (k if k >= 0 else -1 - k)
    assert image >= 0


def check_base():
    monomials = [t**k for k in range(-7, 8)]
    samples = [0, *monomials, (1 + t) / t**2, (1 - t**2) / (1 - t), 3 * t / (2 + t**3)]
    for sample in samples:
        check_kochen(sample, (t,), (1,))
    print(f"Proposition 6.2(1–2): {len(samples)} exact samples in Q(t)")


def check_box(xs):
    # Choose a point avoiding this finite list of zero sets inside (1, 2)^d.
    polynomials = [(3 * x - 4) * (2 * x - 3) for x in xs]
    polynomials.append(sum(xs) - Rational(3 * len(xs), 2))
    grid = tuple(Rational(n, 6) for n in range(7, 12))
    for point in product(grid, repeat=len(xs)):
        substitution = dict(zip(xs, point, strict=True))
        if all(p.subs(substitution) != 0 for p in polynomials):
            assert all(1 < x < 2 for x in point)
            return point
    raise AssertionError("finite box fixture has no witness")


def check_family():
    for d in range(1, len(variables) + 1):
        xs = tuple(variables[:d])
        first = xs[0]
        generators = (t, *xs)
        weights = (1, -1, *((0,) * (d - 1)))
        assert [valuation(x, generators, weights) for x in (1, 2, t, 1 / t, first)] == [
            0, 0, 1, -1, -1
        ]
        samples = [
            0, 1, t, 1 / t, first, t * first, first / t, 1 / (1 + first),
            (1 + t * first) / (2 - t * first),
            (t * first**2 + sum(xs)) / (1 + t**2 * first**2),
            (1 + first) / (1 + t**2 * sum(xs)),
        ]
        for sample in samples:
            check_kochen(sample, generators, weights)
            assert valuation(sample, generators, weights) == valuation(
                cancel(sympify(sample).subs(first, first / t)), generators, (1, *((0,) * d))
            )
        point = check_box(xs)
        print(f"d={d}: {len(samples)} Kochen values, w(1/t)=w(u1)=-1, box witness {point}")


if __name__ == "__main__":
    check_base()
    check_family()
    print("ALL_CHECKS_OK")
