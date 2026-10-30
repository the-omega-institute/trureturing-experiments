#!/usr/bin/env python3
"""Exact finite checks for Joshi's thesis Conjecture 6.1.2.

Qsqrt implements Q(sqrt(D)) using fractions; signs use rational square
comparisons. The independent circle computation intersects open intervals
of all possible initial phases. Irrational-orbit density identifies a
nonempty open phase interval with actual starts in the source word.
No finite prefix is used as an upper bound; no float is used as evidence.
These finite checks do not prove the infinite family or density bridge.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import total_ordering
from math import isqrt
import platform


def require(condition, message):
    if not condition:
        raise AssertionError(message)


@total_ordering
@dataclass(frozen=True, eq=False)
class Qsqrt:
    D: int
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))
        require(self.D > 0 and isqrt(self.D) ** 2 != self.D, 'positive nonsquare radicand')

    def coerce(self, other):
        if isinstance(other, Qsqrt):
            require(self.D == other.D, 'same quadratic field')
            return other
        return Qsqrt(self.D, F(other))

    def __add__(self, other):
        other = self.coerce(other)
        return Qsqrt(self.D, self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Qsqrt(self.D, -self.a, -self.b)

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        return Qsqrt(self.D, self.a * other.a + self.D * self.b * other.b,
                     self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def inverse(self):
        norm = self.a * self.a - self.D * self.b * self.b
        require(norm != 0, 'nonzero denominator')
        return Qsqrt(self.D, self.a / norm, -self.b / norm)

    def __truediv__(self, other):
        return self * self.coerce(other).inverse()

    def __rtruediv__(self, other):
        return self.coerce(other) * self.inverse()

    def __pow__(self, n):
        require(isinstance(n, int) and n >= 0, 'natural exponent')
        result, base = Qsqrt(self.D, 1), self
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n //= 2
        return result

    def sign(self):
        a, b = self.a, self.b
        if b == 0:
            return (a > 0) - (a < 0)
        if a == 0:
            return (b > 0) - (b < 0)
        if a > 0 and b > 0:
            return 1
        if a < 0 and b < 0:
            return -1
        z = a * a - self.D * b * b
        s = (z > 0) - (z < 0)
        return s if a > 0 else -s

    def __eq__(self, other):
        other = self.coerce(other)
        return self.a == other.a and self.b == other.b

    def __lt__(self, other):
        return (self - other).sign() < 0

    def __le__(self, other):
        return (self - other).sign() <= 0

    def __abs__(self):
        return -self if self < 0 else self

    def floor(self):
        k = isqrt(self.D)
        lower = self.a + self.b * (k if self.b >= 0 else k + 1)
        upper = self.a + self.b * (k + 1 if self.b >= 0 else k)
        lo = lower.numerator // lower.denominator - 1
        hi = upper.numerator // upper.denominator + 2
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self < mid:
                hi = mid
            else:
                lo = mid
        require(lo <= self < lo + 1, 'floor certificate')
        return lo

    def ceil(self):
        return -(-self).floor()

    def fract(self):
        return self - self.floor()


def circle_maximum(alpha, d, cap=5000):
    """Intersect I - k*d*alpha modulo 1 until empty, for each symbol I.

    Endpoints cannot be actual source phases, by irrationality. Isolated
    endpoint contacts are discarded. Exceeding the cap raises an error.
    """
    zero, one = Qsqrt(alpha.D), Qsqrt(alpha.D, 1)
    delta = (d * alpha).fract()
    require(0 < delta < 1, 'irrational nonzero rotation step')
    maxima = []
    for L, R in [(zero, one - alpha), (one - alpha, one)]:
        feasible, shift = [(zero, one)], zero
        for length in range(1, cap + 1):
            left, right = L - shift, R - shift
            if right <= 0:
                moved = [(left + 1, right + 1)]
            elif left < 0:
                moved = [(zero, right), (left + 1, one)]
            else:
                moved = [(left, right)]
            intersection = []
            for l, r in feasible:
                for u, v in moved:
                    a, b = max(l, u), min(r, v)
                    if a < b:
                        intersection.append((a, b))
            if not intersection:
                maxima.append(length - 1)
                break
            feasible = intersection
            shift = (shift + delta).fract()
        else:
            raise AssertionError('circle cap reached without a global upper bound')
    return max(maxima), tuple(maxima)


def convergents(digits):
    pm2, pm1, qm2, qm1 = 1, 0, 0, 1
    yield 0, 0, 1
    for n, digit in enumerate(digits, 1):
        p, q = digit * pm1 + pm2, digit * qm1 + qm2
        yield n, p, q
        pm2, pm1, qm2, qm1 = pm1, p, qm1, q


def source_control():
    alpha = Qsqrt(5, -2, 1)
    data = list(convergents([4] * 5))
    require([q for _, _, q in data[:4]] == [1, 4, 17, 72], 'source q indexing')
    for n, p, q in data[1:]:
        require(q * alpha - p == (-1) ** n * alpha ** (n + 1), 'Proposition 4.4.5')
    for d, expected in [(17, 59), (34, 30), (51, 20), (94316, 61)]:
        observed, symbols = circle_maximum(alpha, d)
        require(observed == expected, f'source example A({d})')
        print(f'FIDELITY alpha=sqrt(5)-2 d={d} all_phase_A={observed} letters={symbols}')
    print('FIDELITY q_0,q_1,q_2,q_3=1,4,17,72; constant-period errors checked n=1..5')


def family_checks():
    recurrence_count = circle_count = multiplier_count = 0
    for a in range(3, 9):
        D = a * a + 4 * a
        require((a + 1) ** 2 < D < (a + 2) ** 2, 'irrational family certificate')
        alpha = Qsqrt(D, F(-1, 2), F(1, 2 * a))
        require(0 < alpha < F(1, a) < F(1, 2), 'slope range')
        require(a * alpha * alpha + a * alpha == 1, 'family equation')
        tail = alpha
        for digit in [a, 1] * 4:
            reciprocal = tail.inverse()
            require(reciprocal.floor() == digit, 'literal continued-fraction digits')
            tail = reciprocal.fract()
        require(tail == alpha, 'period-two return')
        for n, p, q in list(convergents([a, 1] * 4))[1:]:
            m = a ** ((n + 1) // 2)
            E = m * alpha ** (n + 1)
            require(q * alpha - p == (-1) ** n * E, 'exact convergent recurrence')
            require(0 < E < alpha, 'small-step condition')
            norm = min((q * alpha).fract(), 1 - (q * alpha).fract())
            require(norm == E, 'nearest-integer error')
            X = (1 - alpha) / E
            actual = X.ceil()
            predicted = ((1 - alpha) / alpha ** (n + 1)).ceil()
            require(X > 1 and (m - 1) * X > 1, 'strict ceiling-gap mechanism')
            require(actual < predicted, 'refutation in checked family member')
            recurrence_count += 1
            if n <= 4:
                observed, symbols = circle_maximum(alpha, q)
                require(observed == actual, 'independent all-phase maximum')
                circle_count += 1
                print(f'FAMILY a={a} n={n} p={p} q={q} all_phase_A={observed}'
                      f' letters={symbols} conjectured={predicted} correction_factor={m}')
                for y in range(2, 5):
                    if y * E <= alpha:
                        require(circle_maximum(alpha, y * q)[0] == ((1 - alpha) / (y * E)).ceil(),
                                'multiplier correction')
                        multiplier_count += 1
    alpha = Qsqrt(21, F(-1, 2), F(1, 6))
    require((1 - alpha) / (3 * alpha ** 2) == Qsqrt(21, 2, F(1, 3)), 'actual argument')
    require((1 - alpha) / alpha ** 2 == Qsqrt(21, 6, 1), 'conjectured argument')
    letters = [1 - alpha <= ((i + 1) * alpha).fract() for i in [13, 16, 19, 22]]
    require(letters == [False] * 4, 'source-indexed four-letter occurrence')
    print('WITNESS a=3,n=1: ceil(2+sqrt(21)/3)=4 < ceil(6+sqrt(21))=11;'
          ' source positions 13,16,19,22 all carry the first letter')
    print(f'COUNTS exact_family_recurrences={recurrence_count}'
          f' independent_family_circle_maxima={circle_count}'
          f' admissible_multiplier_circle_maxima={multiplier_count}')
    print('SCOPE finite exact checks; the universal family and density bridge require proof')


def main():
    print(f'Python {platform.python_version()}; integers and fractions; no floating-point evidence')
    source_control()
    family_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
