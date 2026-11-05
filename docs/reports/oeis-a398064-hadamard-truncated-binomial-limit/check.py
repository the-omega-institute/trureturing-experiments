#!/usr/bin/env python3
"""Exact checks for OEIS A398064; Python 3.10+ and SymPy 1.14.

No binary floating point or decimal approximation is used. The finite checks and
symbolic limits below are evidence for the preregistration, not a proof over C.
"""
from fractions import Fraction
from math import comb, factorial
import sympy as sp


def signed_binomial(x, k):
    """(-1)^k x(x-1)...(x-k+1)/k!, computed by rational recurrence."""
    c = Fraction(1)
    for j in range(k):
        c *= (j - x) / (j + 1)
    return c


def rational_sum(x, n):
    assert n >= 1
    return sum((signed_binomial(x, k) / (n + k) for k in range(n + 1)), Fraction())


def rational_scaled_integer(x, n):
    return Fraction(n) ** (x + 1) * rational_sum(Fraction(x), n)


def exact_zero(expr):
    assert not expr.has(sp.Float)
    assert sp.simplify(expr) == 0, expr


def main():
    # Source fidelity: its two exact formulas at the quarter argument agree.
    q = sp.Rational(1, 4)
    digamma_form = (sp.digamma(sp.Rational(7, 8)) - sp.digamma(sp.Rational(3, 8))) / (2 * sp.gamma(1 - q))
    radical_form = (sp.pi + sp.log(3 - 2 * sp.sqrt(2))) / (sp.sqrt(2) * sp.gamma(sp.Rational(3, 4)))
    exact_zero(sp.expand_log(sp.expand_func(digamma_form) - radical_form, force=True))
    assert sp.simplify((2 - sp.sqrt(2)) / (2 + sp.sqrt(2)) - (3 - 2 * sp.sqrt(2))) == 0
    assert 3 - 2 * sp.sqrt(2) > 0
    print('FIDELITY: quarter-value digamma/radical identity verified exactly')

    # All nonnegative integers in the declared small-case window.
    N = sp.Symbol('N', integer=True, positive=True)
    nonnegative_cases = 0
    for m in range(9):
        finite = N ** (m + 1) * sum(sp.Rational((-1) ** k * comb(m, k), 1) / (N + k) for k in range(m + 1))
        beta_form = sp.factorial(m) * N ** (m + 1) / sp.prod(N + j for j in range(m + 1))
        exact_zero(sp.cancel(finite - beta_form))
        assert sp.limit(beta_form, N, sp.oo) == sp.factorial(m)
        for n in range(max(1, m), 65):
            lhs = rational_scaled_integer(m, n)
            rhs = Fraction(factorial(m) * n ** (m + 1), 1)
            for j in range(m + 1):
                rhs /= n + j
            assert lhs == rhs
            nonnegative_cases += 1
    print(f'INTEGER LIMITS: x=0,...,8 give x! exactly; {nonnegative_cases} finite-sum checks')

    # For x=-p the coefficients are polynomials. Divide by N+k, sum the
    # polynomial part exactly, and keep h=H_(2N)-H_(N-1) as a symbolic variable.
    # h -> log(2); both coefficients have finite limits, checked separately.
    k, h = sp.symbols('k h')
    beta_integer = sp.log(2)
    negative_cases = 0
    for p in range(1, 9):
        polynomial = sp.prod(k + j for j in range(1, p)) / sp.factorial(p - 1)
        quotient, remainder = sp.div(polynomial, k + N, k)
        poly_sum = sum(coef * sp.summation(k ** power[0], (k, 0, N))
                       for power, coef in sp.Poly(quotient, k).terms())
        closed = sp.expand(N ** (1 - p) * (poly_sum + remainder * h))
        a = sp.expand(closed).coeff(h, 1)
        b = sp.expand(closed).subs(h, 0)
        a_limit = sp.limit(a, N, sp.oo)
        b_limit = sp.limit(b, N, sp.oo)
        assert a_limit.is_finite and b_limit.is_finite
        exact_zero(a_limit * sp.log(2) + b_limit - beta_integer / sp.factorial(p - 1))
        for n in range(1, 65):
            hn = sum((Fraction(1, n + j) for j in range(n + 1)), Fraction())
            symbolic = closed.subs({N: n, h: sp.Rational(hn.numerator, hn.denominator)})
            lhs = rational_scaled_integer(-p, n)
            assert symbolic == sp.Rational(lhs.numerator, lhs.denominator)
            negative_cases += 1
        beta_integer = sp.Rational(1, p) - beta_integer
    print(f'INTEGER LIMITS: x=-1,...,-8 give beta(-x)/Gamma(-x) exactly; {negative_cases} finite-sum checks')

    # Rational interval for log(2)=2 atanh(1/3). The remaining positive series
    # is bounded by a geometric tail with denominators at least 2P+1.
    P = 32
    log_lo = sum((Fraction(2, (2 * j + 1) * 3 ** (2 * j + 1)) for j in range(P)), Fraction())
    log_hi = log_lo + Fraction(2, (2 * P + 1) * 3 ** (2 * P + 1)) / (1 - Fraction(1, 9))
    for n in range(1, 129):
        s1 = rational_scaled_integer(-1, n)
        assert Fraction(1, 2 * n) <= s1 - log_hi
        assert s1 - log_lo <= Fraction(1, n)
        s2 = rational_scaled_integer(-2, n)
        assert -Fraction(2, n) <= s2 - (1 - log_lo)
        assert s2 - (1 - log_hi) <= Fraction(2, n)
    print('RATIONAL INTERVALS: x=-1 and x=-2 satisfy explicit O(1/n) bounds for n=1,...,128')

    # Exact polynomial identity behind the Euler GammaSeq coefficient scaling.
    x = sp.Symbol('x')
    for j in range(13):
        cj = (-1) ** j * sp.prod(x - r for r in range(j)) / sp.factorial(j)
        cj1 = (-1) ** (j + 1) * sp.prod(x - r for r in range(j + 1)) / sp.factorial(j + 1)
        exact_zero((j + 1) * cj1 - (j - x) * cj)
        exact_zero(sp.prod(r - x for r in range(j + 1)) - (j - x) * sp.factorial(j) * cj)
    print('COEFFICIENTS: recurrence and GammaSeq product identity verified symbolically for k=0,...,12')

    # Separate product and recurrence implementations on rational and genuinely
    # complex inputs; these check algebra and conventions, not analytic limits.
    values = [sp.Rational(-7, 4), sp.Rational(-3, 4), sp.Rational(-1, 2),
              sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(3, 2),
              (1 + sp.I) / 2, -1 + sp.I]
    algebra_cases = 0
    for z in values:
        coeffs = [sp.Integer(1)]
        for j in range(12):
            coeffs.append(sp.expand(coeffs[-1] * (j - z) / (j + 1)))
        for n in range(1, 13):
            lhs = sum(coeffs[j] / (n + j) for j in range(n + 1))
            rhs = sum((-1) ** j * sp.binomial(z, j) / (n + j) for j in range(n + 1))
            exact_zero(sp.expand_func(lhs - rhs))
            algebra_cases += 1
    print(f'RATIONAL/COMPLEX CONVENTIONS: {algebra_cases} exact finite-sum cross-checks')
    print('SCOPE: symbolic limits cover integer x=-8,...,8 only; no full complex-plane limit is certified')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
