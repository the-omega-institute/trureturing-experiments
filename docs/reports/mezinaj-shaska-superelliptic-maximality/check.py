#!/usr/bin/env python3
"""Exact checks for Mezinaj–Shaska Problem 7.11; Python standard library only.
Field arithmetic uses integer coefficient vectors modulo a verified irreducible
polynomial. No floating point, probabilistic primality test, or external data.
"""
from itertools import product
from collections import Counter
from fractions import Fraction
from math import gcd


def trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def rem(a, b, p):
    a = trim(list(a))
    while len(a) >= len(b):
        c = a[-1] * pow(b[-1], -1, p) % p
        j = len(a) - len(b)
        for i, v in enumerate(b):
            a[i + j] = (a[i + j] - c * v) % p
        trim(a)
    return a


def irreducible(f, p):
    # A reducible degree-m polynomial has a monic factor of degree <= m/2.
    m = len(f) - 1
    return all(rem(f, list(c) + [1], p)
               for d in range(1, m // 2 + 1)
               for c in product(range(p), repeat=d))


class Field:
    def __init__(self, p, m):
        assert p >= 2 and all(p % k for k in range(2, p))
        self.p, self.m, self.q = p, m, p ** m
        self.modulus = next(list(c) + [1] for c in product(range(p), repeat=m)
                            if c[0] and irreducible(list(c) + [1], p))
        self.vec = [tuple((a // (p ** i)) % p for i in range(m))
                    for a in range(self.q)]

    def encode(self, v):
        return sum((a % self.p) * self.p ** i for i, a in enumerate(v))

    def add(self, a, b):
        return self.encode([x + y for x, y in zip(self.vec[a], self.vec[b])])

    def neg(self, a):
        return self.encode([-x for x in self.vec[a]])

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        c = [0] * (2 * self.m - 1)
        for i, x in enumerate(self.vec[a]):
            for j, y in enumerate(self.vec[b]):
                c[i + j] += x * y
        return self.encode(rem(c, self.modulus, self.p))

    def power(self, a, e):
        assert e >= 0
        b = 1
        while e:
            if e % 2:
                b = self.mul(a, b)
            a = self.mul(a, a)
            e //= 2
        return b

    def eval(self, f, a):
        out = 0
        for c in reversed(f):
            out = self.add(self.mul(out, a), c)
        return out

    def derivative(self, f):
        return [self.mul(i % self.p, f[i]) for i in range(1, len(f))]

    def polynomial_product(self, f, g):
        out = [0] * (len(f) + len(g) - 1)
        for i, x in enumerate(f):
            for j, y in enumerate(g):
                out[i + j] = self.add(out[i + j], self.mul(x, y))
        return trim(out)

    def verify(self):
        # Independent exhaustive inverse, Frobenius and quadratic-fiber checks.
        for a in range(self.q):
            assert self.power(a, self.q) == a
            if a:
                assert self.mul(a, self.power(a, self.q - 2)) == 1
        squares = Counter(self.mul(y, y) for y in range(self.q))
        assert squares[0] == 1
        assert all(v == 2 for a, v in squares.items() if a)


def genus(n, d):
    g = Fraction((n - 1) * (d - 1) + 1 - gcd(n, d), 2)
    assert g.denominator == 1 and g >= 0
    return g.numerator


def fidelity(F, r):
    # Source Theorem 7.10, including all n >= 2 dividing r+1 in these fields.
    for n in range(2, r + 2):
        if (r + 1) % n:
            continue
        values = [F.neg(F.add(F.power(a, r + 1), 1)) for a in range(F.q)]
        fibers = Counter(F.power(y, n) for y in range(F.q))
        roots = sum(v == 0 for v in values)
        nonbranch = sum(fibers[v] for v in values if v)
        affine = sum(fibers[v] for v in values)
        infinity = sum(F.power(y, gcd(n, r + 1)) == F.neg(1)
                       for y in range(F.q))
        delta = sum(v != 0 and fibers[v] == 0 for v in values)
        g = genus(n, r + 1)
        assert roots == r + 1 and delta == 0
        assert nonbranch == n * (r * r - r - 1)
        assert infinity == n
        assert g == (n - 1) * (r - 1) // 2
        assert affine + infinity == F.q + 1 + 2 * r * g
        print(f"FIDELITY r={r} n={n} d={roots} delta={delta} "
              f"N_evaluation={nonbranch} affine={affine} infinity={infinity} "
              f"genus={g} projective={affine + infinity} maximal=True")


def counterexample(F, r):
    q, d = r * r, r * r - r
    assert q == F.q and r >= 3 and r % 2 == 1
    assert F.p % 2 == 1 and (q - 1) % 2 == 0 and d >= 6 and d % 2 == 0
    f = [0] * (d + 1)
    for j in range(r + 1):
        f[j * (r - 1)] = 1
    # Verify the polynomial quotient identity without division or value cancellation.
    denominator = [0] * (r + 1)
    denominator[1], denominator[r] = F.neg(1), 1
    numerator = [0] * (q + 1)
    numerator[1], numerator[q] = F.neg(1), 1
    assert F.polynomial_product(denominator, f) == numerator
    subfield = [a for a in range(q) if F.power(a, r) == a]
    assert len(subfield) == r
    complement = [a for a in range(q) if a not in subfield]
    # A second construction, as a product of all complementary linear factors.
    root_product = [1]
    for a in complement:
        root_product = F.polynomial_product(root_product, [F.neg(a), 1])
    assert root_product == f
    values = [F.eval(f, a) for a in range(q)]
    assert all(v == (1 if a in subfield else 0) for a, v in enumerate(values))
    deriv = F.derivative(f)
    assert all(F.eval(deriv, a) != 0 for a in complement)
    powers = Counter(F.power(y, 2) for y in range(q))
    affine_points = [(a, y) for a, v in enumerate(values)
                     for y in range(q) if F.power(y, 2) == v]
    affine = len(affine_points)
    assert affine == sum(powers[v] for v in values) == q + r
    delta = sum(v != 0 and powers[v] == 0 for v in values)
    assert delta == 0
    # Reverse-polynomial chart v^2=t^d f(1/t); no inverse is evaluated at zero.
    reverse = list(reversed(f))
    infinity_points = [v for v in range(q) if F.power(v, 2) == F.eval(reverse, 0)]
    assert len(infinity_points) == 2
    assert all(F.mul(2 % F.p, v) != 0 for v in infinity_points)
    for a, y in affine_points:
        if not a:
            continue
        t = F.power(a, q - 2)
        v = F.mul(y, F.power(t, d // 2))
        assert F.power(v, 2) == F.eval(reverse, t)
        assert F.mul(t, a) == 1 and F.mul(v, F.power(a, d // 2)) == y
    projective = affine + len(infinity_points)
    g = genus(2, d)
    bound = q + 1 + 2 * r * g
    assert projective == q + r + 2
    assert 0 < r + 1 < 2 * r
    assert (projective - q - 1) % (2 * r) != 0
    assert bound - projective == (r + 1) * (r * r - 2 * r - 1) > 0
    print(f"COUNTEREXAMPLE p={F.p} r={r} q={q} modulus={F.modulus} "
          f"d={d} roots={len(complement)} nonroots={len(subfield)} "
          f"delta={delta} affine={affine} infinity={len(infinity_points)} "
          f"projective={projective} genus={g} HW_maximum={bound} gap={bound-projective}")
    if r == 3:
        print("F9 coordinate vector [a,b] means a+b*i, i^2=2")
        print("F9 evaluations=" + repr([(F.vec[a], F.vec[v]) for a, v in enumerate(values)]))
        print("F9 root_derivatives=" + repr([(F.vec[a], F.vec[F.eval(deriv,a)])
                                             for a in complement]))
        print("F9 affine_points=" + repr([(F.vec[a], F.vec[y]) for a, y in affine_points]))
        print("F9 infinity=" + repr([F.vec[v] for v in infinity_points]))


def main():
    for p, a in [(3, 1), (5, 1), (7, 1), (3, 2), (11, 1), (13, 1)]:
        r = p ** a
        F = Field(p, 2 * a)
        F.verify()
        fidelity(F, r)
        counterexample(F, r)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
