#!/usr/bin/env python3
"""Exact finite checks for Ferreira--Simon, arXiv:2108.09244v2, Section 5.5.

Every interval endpoint is an integer divided by SCALE. Arithmetic rounds
outwards. No floating-point operations or external packages are used.
These finite checks are not a proof of the all-index analytic argument.
"""
from fractions import Fraction as Q
from math import factorial, isqrt

SCALE = 10 ** 240


def ceildiv(a, b):
    assert b > 0
    return -((-a) // b)


class I:
    """Closed interval [lo/SCALE, hi/SCALE]."""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi):
        assert lo <= hi
        self.lo, self.hi = lo, hi

    @staticmethod
    def rational(x):
        x = Q(x)
        return I((x.numerator * SCALE) // x.denominator,
                 ceildiv(x.numerator * SCALE, x.denominator))

    @staticmethod
    def coerce(x):
        return x if isinstance(x, I) else I.rational(x)

    def __add__(self, other):
        other = I.coerce(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -I.coerce(other)

    def __rsub__(self, other):
        return I.coerce(other) + -self

    def __mul__(self, other):
        other = I.coerce(other)
        p = [a * b for a in (self.lo, self.hi)
             for b in (other.lo, other.hi)]
        return I(min(p) // SCALE, ceildiv(max(p), SCALE))

    __rmul__ = __mul__

    def reciprocal(self):
        assert self.lo > 0 or self.hi < 0, 'division interval meets zero'
        # a/b with negative b is reduced to -a/(-b).
        def floor_ratio(a, b):
            return a // b
        def ceil_ratio(a, b):
            return -((-a) // b)
        return I(floor_ratio(SCALE * SCALE, self.hi),
                 ceil_ratio(SCALE * SCALE, self.lo))

    def __truediv__(self, other):
        if not isinstance(other, I):
            q = Q(other)
            assert q != 0
            a, b = self.lo * q.denominator, self.hi * q.denominator
            if q < 0:
                a, b = -b, -a
            d = abs(q.numerator)
            return I(a // d, ceildiv(b, d))
        return self * other.reciprocal()

    def __rtruediv__(self, other):
        return I.coerce(other) / self

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        result, base = I.rational(1), self
        while n:
            if n % 2:
                result = result * base
            base = base * base
            n //= 2
        return result

    def inside(self, lo, hi):
        return Q(lo) < Q(self.lo, SCALE) and Q(self.hi, SCALE) < Q(hi)


def atan_bounds(d, terms=220):
    """Alternating arctan(1/d) series and its first omitted term."""
    s = sum((Q((-1) ** j, (2 * j + 1) * d ** (2 * j + 1))
             for j in range(terms)), Q(0))
    t = Q((-1) ** terms, (2 * terms + 1) * d ** (2 * terms + 1))
    return min(s, s + t), max(s, s + t)


def pi_bounds():
    a, b = atan_bounds(5), atan_bounds(239)
    # Machin's identity: pi = 16 atan(1/5) - 4 atan(1/239).
    return 16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0]


def sqrt_bounds(lo, hi):
    assert 0 <= lo <= hi
    l = isqrt((lo.numerator * SCALE * SCALE) // lo.denominator)
    h = isqrt((hi.numerator * SCALE * SCALE) // hi.denominator) + 1
    assert Q(l * l, SCALE * SCALE) <= lo
    assert hi <= Q(h * h, SCALE * SCALE)
    return I(l, h)


def exp_bounds(y, terms=260):
    assert 0 <= y <= Q(1, 2)
    s = sum((y ** j / factorial(j) for j in range(terms)), Q(0))
    first = y ** terms / factorial(terms)
    # Ratios of successive omitted terms are at most y/(terms+1).
    tail = first / (1 - y / (terms + 1))
    return I(I.rational(s).lo, I.rational(s + tail).hi)


def mills(a):
    """Equation defining r: exp(a²/2) [sqrt(pi/2) - int_0^a exp(-t²/2)dt]."""
    assert 0 <= a <= 1
    terms = 260
    s = sum((Q((-1) ** j) * a ** (2 * j + 1) /
             (2 ** j * factorial(j) * (2 * j + 1))
             for j in range(terms)), Q(0))
    t = Q((-1) ** terms) * a ** (2 * terms + 1) / (
        2 ** terms * factorial(terms) * (2 * terms + 1))
    integral = I(I.rational(min(s, s + t)).lo,
                 I.rational(max(s, s + t)).hi)
    lo, hi = pi_bounds()
    return exp_bounds(a * a / 2) * (sqrt_bounds(lo / 2, hi / 2) - integral)


def derivatives(a, max_order):
    """Differentiated Mills ODE: D[j+1] = a D[j] + j D[j-1]."""
    d = [mills(a)]
    if max_order:
        d.append(a * d[0] - 1)
    for j in range(1, max_order):
        d.append(a * d[j] + j * d[j - 1])
    # The Laplace integral gives (-1)^j D[j] > 0 on the real axis.
    assert all(((-1) ** j * v).lo > 0 for j, v in enumerate(d))
    return d


def convolution(a, b, degree):
    return [sum((a[j] * b[k - j] for j in range(k + 1)), I.rational(0))
            for k in range(degree + 1)]


def quotient_jet(n, a, degree):
    d = derivatives(a, 2 * n + 1 + degree)
    u = [d[n + k] / factorial(k) for k in range(degree + 1)]
    b = [d[2 * n + 1 + k] / factorial(k) for k in range(degree + 1)]
    numerator = [-v for v in convolution(u, u, degree)]
    q = []
    for k in range(degree + 1):
        q.append((numerator[k] - sum(
            (b[j] * q[k - j] for j in range(1, k + 1)), I.rational(0))) / b[0])
    return q


def sqrt_composition(q):
    """Jet of F(sqrt(z)) at z=1 from the jet of F at x=1."""
    degree = len(q) - 1
    # sqrt(1+h)-1 = sum binom(1/2,k) h^k for k >= 1.
    u = [I.rational(0)]
    c = Q(1)
    for k in range(1, degree + 1):
        c *= (Q(1, 2) - k + 1) / k
        u.append(I.rational(c))
    power = [I.rational(int(k == 0)) for k in range(degree + 1)]
    result = [I.rational(0) for _ in range(degree + 1)]
    for coefficient in q:
        result = [x + coefficient * y for x, y in zip(result, power)]
        power = convolution(power, u, degree)
    return result


def compact_negative_bound(v):
    """A short rational enclosure of a strictly negative interval."""
    assert v.hi < 0
    magnitude = 0
    bound = Q(-v.hi, SCALE)
    while bound < Q(1, 10 ** magnitude):
        magnitude += 1
    digits = magnitude + 5
    scale = 10 ** digits
    l = (v.lo * scale) // SCALE - 1
    h = ceildiv(v.hi * scale, SCALE) + 1
    assert v.inside(Q(l, scale), Q(h, scale)) and h < 0
    return f'{l}/10^{digits} < coefficient < {h}/10^{digits} < 0'


def main():
    print('Arithmetic: integers and fractions; outward scale = 10^240; no floats')
    # Source fidelity: Sampford corollary in Section 5.5, at x=0, 1/100, 1.
    for a in (Q(0), Q(1, 100), Q(1)):
        s = sqrt_bounds(a * a + 8, a * a + 8)
        gap = 4 / (3 * a + s) - mills(a)
        assert gap.lo > 0
    print('FIDELITY: r(x) < 4/(3x+sqrt(x^2+8)), x = 0, 1/100, 1')
    # A finite derivative check of Proposition 10, NOT a proof of complete monotonicity.
    for n in range(5):
        composed = sqrt_composition(quotient_jet(n, Q(1), 12))
        assert all(((-1) ** k * v).lo > 0 for k, v in enumerate(composed))
    print('FIDELITY: Proposition 10 alternating derivatives at z=1, n=0..4, k=0..12')
    a = Q(1, 100)
    jets = [quotient_jet(n, a, 220) for n in range(5)]
    # The two independently supplied witnesses are reproduced with strict rational bounds.
    assert jets[0][40].inside(Q(-2995, 10 ** 29), Q(-2994, 10 ** 29))
    assert jets[1][70].inside(Q(-4604, 10 ** 55), Q(-4603, 10 ** 55))
    print('F_0^(40)(1/100)/40!: -2995/10^29 < coefficient < -2994/10^29 < 0')
    print('F_1^(70)(1/100)/70!: -4604/10^55 < coefficient < -4603/10^55 < 0')
    for n, jet in enumerate(jets):
        negatives = [(k, (-1) ** k * v) for k, v in enumerate(jet)
                     if ((-1) ** k * v).hi < 0]
        assert negatives, ('no negative alternating coefficient found', n)
        k, v = negatives[0]
        print(f'REFUTATION SAMPLE n={n}, x=1/100, k={k}: ' + compact_negative_bound(v))
    print('Scope: finite exact witnesses for n=0..4; all-index/every-center argument is analytic')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
