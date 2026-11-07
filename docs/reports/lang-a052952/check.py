#!/usr/bin/env python3
"""Exact checks for the complex-Fibonacci floor conjecture in OEIS A052952.

Only integers and fractions.Fraction enter the calculations. Each analytic
quantity is enclosed by rational bounds with explicit series remainder bounds.
The finite check is evidence for the preregistered route, not an all-index proof.
"""

from dataclasses import dataclass
from fractions import Fraction as Q
from math import isqrt


SCALE = 1 << 384


def down(x):
    return Q((x * SCALE).__floor__(), SCALE)


def up(x):
    return Q((x * SCALE).__ceil__(), SCALE)


@dataclass(frozen=True)
class I:
    lo: Q
    hi: Q

    def __post_init__(self):
        assert self.lo <= self.hi

    @staticmethod
    def point(x):
        return I(Q(x), Q(x))

    @staticmethod
    def enclose(lo, hi):
        return I(down(Q(lo)), up(Q(hi)))

    def __add__(self, other):
        other = interval(other)
        return I.enclose(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -interval(other)

    def __rsub__(self, other):
        return interval(other) - self

    def __mul__(self, other):
        other = interval(other)
        products = [a * b for a in (self.lo, self.hi)
                    for b in (other.lo, other.hi)]
        return I.enclose(min(products), max(products))

    __rmul__ = __mul__

    def inverse(self):
        assert self.lo > 0 or self.hi < 0
        return I.enclose(1 / self.hi, 1 / self.lo)

    def __truediv__(self, other):
        return self * interval(other).inverse()

    def __rtruediv__(self, other):
        return interval(other) / self

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        ans, base = I.point(1), self
        while n:
            if n % 2:
                ans = ans * base
            base = base * base
            n //= 2
        return ans


def interval(x):
    return x if isinstance(x, I) else I.point(x)


def arctan_reciprocal(q, terms=180):
    """Alternating-series bracket for arctan(1/q), q > 1."""
    assert q > 1
    s = sum((Q((-1) ** k, (2 * k + 1) * q ** (2 * k + 1))
             for k in range(terms)), Q(0))
    next_term = Q((-1) ** terms, (2 * terms + 1) * q ** (2 * terms + 1))
    return I.enclose(min(s, s + next_term), max(s, s + next_term))


def log_positive(x, terms=160):
    """log x = 2 sum z^(2k+1)/(2k+1), z=(x-1)/(x+1)."""
    assert x.lo > 1
    z = (x - 1) / (x + 1)
    assert 0 < z.lo <= z.hi < 1
    term, s = z, I.point(0)
    for k in range(terms):
        s = s + term / (2 * k + 1)
        term = term * z * z
    # All omitted terms are positive. Bound their denominators from below
    # by 2*terms+1 and sum the remaining geometric series.
    tail = 2 * term / ((2 * terms + 1) * (1 - z * z))
    value = 2 * s
    return I.enclose(value.lo, value.hi + tail.hi)


def exp_small_positive(x, terms=100):
    """Positive exponential series, with a geometric bound on its tail."""
    assert 0 <= x.lo <= x.hi < 1
    term, s = I.point(1), I.point(0)
    for k in range(terms):
        s = s + term
        term = term * x / (k + 1)
    # term encloses x^terms/terms!; later successive ratios are at most
    # x/(terms+1). The whole omitted tail includes this first omitted term.
    tail = term / (1 - x / (terms + 1))
    return I.enclose(s.lo, s.hi + tail.hi)


def cosine(x, terms=180):
    """Taylor polynomial through degree 2*terms-1, with Lagrange remainder."""
    term, s = I.point(1), I.point(0)
    for k in range(terms):
        s = s + term
        term = -term * x * x / ((2 * k + 1) * (2 * k + 2))
    # |cos^(2*terms)| <= 1 on the real line. The next even term encloses
    # x^(2*terms)/(2*terms)! up to its sign, including interval roundoff.
    error = max(abs(term.lo), abs(term.hi))
    return s + I(-error, error)


def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def main():
    # Fidelity: the 39 displayed terms in revision 166, %S/%T/%U, and
    # entries 39 and 40 in https://oeis.org/A052952/b052952.txt.
    source_terms = [
        1, 1, 3, 4, 8, 12, 21, 33, 55, 88, 144, 232, 377, 609, 987,
        1596, 2584, 4180, 6765, 10945, 17711, 28656, 46368, 75024,
        121393, 196417, 317811, 514228, 832040, 1346268, 2178309,
        3524577, 5702887, 9227464, 14930352, 24157816, 39088169,
        63245985, 102334155, 165580140, 267914296,
    ]
    # The internal entry displays 39 terms through 102334155; the next two
    # are checked against its linked b-file as well as its NAME formula.
    expected = [fib(n + 2) - (1 - (-1) ** n) // 2 for n in range(41)]
    assert source_terms == expected
    # Independently generate coefficients of 1/((1-x-x^2)*(1-x^2)).
    gf_terms = [sum(fib(n - 2 * j + 1) for j in range(n // 2 + 1))
                for n in range(41)]
    assert gf_terms == expected
    for n in range(3, 41):
        assert expected[n] == 2 * expected[n - 2] + expected[n - 3] + 1
    print('FIDELITY_OK: NAME, generating function, recurrence, n=0..40')

    root = isqrt(5 * SCALE * SCALE)
    assert root * root < 5 * SCALE * SCALE < (root + 1) * (root + 1)
    d = I(Q(root, SCALE), Q(root + 1, SCALE))
    phi = (1 + d) / 2
    assert phi.lo > 1 and d.lo > 2
    # Machin's exact identity pi = 16 arctan(1/5) - 4 arctan(1/239).
    pi = 16 * arctan_reciprocal(5) - 4 * arctan_reciprocal(239)
    assert Q(3) < pi.lo < pi.hi < Q(22, 7)
    log_phi = log_positive(phi)
    exp_pi = exp_small_positive(pi / 16) ** 16
    exp_minus_pi = 1 / exp_pi
    assert 0 < exp_minus_pi.lo < exp_minus_pi.hi < Q(1, 2)

    print('RATIONAL_INTERVALS: outward rounding at 2^-384; no float arithmetic')
    print('SERIES: arctan 180; log 160; exp(pi/16) 100; cos 180 terms')
    for n in range(41):
        N, answer = n + 2, expected[n]
        power = phi ** N
        x, delta = power / d, 1 / (power * d)
        q = exp_minus_pi ** N
        c = cosine(2 * N * log_phi)
        # Direct expansion of Lang's definition at N+iN:
        # |F|^2 = x^2 + q^2 delta^2 - 2 (-1)^N q cos(2N log phi)/5.
        # Here x*delta = 1/5 exactly, independently of Binet's formula.
        norm_squared = x * x + q * q * delta * delta - 2 * (-1) ** N * q * c / 5
        assert answer >= 0
        assert answer * answer <= norm_squared.lo
        assert norm_squared.hi < (answer + 1) ** 2
        # A compact rational enclosure of the squared displacement is data,
        # while the assertions above use the full 384-bit enclosures.
        margin = norm_squared - answer * answer
        scale64 = 1 << 64
        lo64 = (margin.lo * scale64).__floor__()
        hi64 = (margin.hi * scale64).__ceil__()
        print(f'n={n:2d} floor={answer} squared_margin=[{lo64}/{scale64},{hi64}/{scale64}]')
    print('FLOORS_0_TO_40=' + ','.join(map(str, expected)))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
