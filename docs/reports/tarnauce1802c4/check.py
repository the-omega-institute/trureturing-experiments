#!/usr/bin/env python3
"""Exact checks of Tarnauceanu C4's proposed coprime-product refutation.

Standard library only. Small group lattices and normal sections are enumerated
from multiplication, independently of the scalar-group counting formula.
Large witnesses and finite samples of the family use integer certificates.
"""
from collections import deque
from fractions import Fraction as F
from itertools import product
from math import gcd, isqrt


class Group:
    def __init__(self, elements, mul):
        self.elements = tuple(elements)
        index = {x: i for i, x in enumerate(self.elements)}
        self.t = [[index[mul(x, y)] for y in self.elements] for x in self.elements]
        self.n = len(self.elements)
        assert all(self.t[0][i] == i == self.t[i][0] for i in range(self.n))
        self.inv = [next(j for j in range(self.n) if self.t[i][j] == 0 == self.t[j][i])
                    for i in range(self.n)]
        self._lattice = None

    def closure(self, generators):
        seen, todo = {0}, [0]
        for x in todo:
            for y in generators:
                z = self.t[x][y]
                if z not in seen:
                    seen.add(z)
                    todo.append(z)
        return frozenset(seen)

    def lattice(self):
        if self._lattice is None:
            known = {frozenset({0})}
            todo = deque(known)
            while todo:
                h = todo.popleft()
                for x in range(self.n):
                    if x not in h:
                        k = self.closure(tuple(h) + (x,))
                        if k not in known:
                            known.add(k)
                            todo.append(k)
            self._lattice = tuple(sorted(known, key=lambda h: (len(h), sorted(h))))
        return self._lattice

    def normal(self, n, h):
        return all(self.t[self.t[x][y]][self.inv[x]] in n for x in h for y in n)

    def sd(self):
        l = self.lattice()
        count = sum({self.t[x][y] for x in h for y in k}
                    == {self.t[y][x] for x in h for y in k} for h in l for k in l)
        return F(count, len(l) ** 2)

    def quotient(self, h, n):
        cosets, owner = [], {}
        for x in sorted(h):
            if x not in owner:
                c = frozenset(self.t[x][y] for y in n)
                for y in c:
                    owner[y] = len(cosets)
                cosets.append(c)
        reps = [min(c) for c in cosets]
        return Group(range(len(cosets)), lambda i, j: owner[self.t[reps[i]][reps[j]]])

    def section_values(self):
        l = self.lattice()
        return {self.quotient(h, n).sd() for h in l for n in l
                if n <= h and self.normal(n, h)}

    def center(self):
        return {x for x in range(self.n)
                if all(self.t[x][y] == self.t[y][x] for y in range(self.n))}

    def proper_subgroups_abelian(self):
        return all(all(self.t[x][y] == self.t[y][x] for x in h for y in h)
                   for h in self.lattice() if len(h) < self.n)


def gaussian_row(p, r):
    # Recurrence, independent of the product formula below.
    row = [1]
    for n in range(1, r + 1):
        row = [1] + [row[j] + p ** (n - j) * row[j - 1]
                     for j in range(1, n)] + [1]
    for j, value in enumerate(row):
        numer, denom = 1, 1
        for i in range(j):
            numer *= p ** (r - i) - 1
            denom *= p ** (j - i) - 1
        assert numer % denom == 0 and value == numer // denom
    return row


def scalar_sd(p, r):
    a, m = sum(gaussian_row(p, r)), p ** r
    return F((a + 1) ** 2 + 7 * m, (a + m + 1) ** 2)


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def order(p, q):
    assert prime(p) and prime(q) and p != q
    x = 1
    for r in range(1, q):
        x = x * p % q
        if x == 1:
            return r
    raise AssertionError('Fermat order missing')


def binary_affine(r, polynomial, q):
    m = 2 ** r
    def mul(x, y):
        z = 0
        while y:
            if y & 1:
                z ^= x
            y >>= 1
            x <<= 1
            if x & m:
                x ^= polynomial
        return z
    # Certify the quotient polynomial defines a field, by all nonzero inverses.
    assert all(any(mul(x, y) == 1 for y in range(1, m)) for x in range(1, m))
    zeta = next(x for x in range(2, m) if field_order(x, mul) == q)
    powers = [1]
    for _ in range(1, q):
        powers.append(mul(powers[-1], zeta))
    return Group(product(range(m), range(q)),
                 lambda x, y: (x[0] ^ mul(powers[x[1]], y[0]), (x[1] + y[1]) % q))


def field_order(x, mul):
    y, n = x, 1
    while y != 1:
        y, n = mul(y, x), n + 1
        assert n < 10000
    return n


def main():
    # Source Theorem 1's exact D8 value, directly from subgroup set products.
    d8 = Group(product(range(4), range(2)),
               lambda x, y: ((x[0] + (-1) ** x[1] * y[0]) % 4, (x[1] + y[1]) % 2))
    assert len(d8.lattice()) == 10 and d8.sd() == F(23, 25)
    assert min(d8.section_values()) == F(23, 25)
    print('SOURCE_FIDELITY D8: subgroups=10 sd=sdStar=23/25')

    a4 = None
    for r, poly, q in [(2, 0b111, 3), (3, 0b1011, 7), (4, 0b10011, 5)]:
        g = binary_affine(r, poly, q)
        assert order(2, q) == r
        a, m = sum(gaussian_row(2, r)), 2 ** r
        assert len(g.lattice()) == a + m + 1
        assert g.sd() == scalar_sd(2, r)
        assert g.section_values() == {F(1), g.sd()}
        assert g.center() == {0} and g.proper_subgroups_abelian()
        assert sum(g.normal(h, frozenset(range(g.n))) for h in g.lattice()) == 3
        if r == 2:
            a4 = g
            assert g.sd() == F(16, 25)  # Source Theorem 2, A4.
        print(f'ENUMERATED_SCALAR p=2 r={r} q={q} order={g.n} '
              f'subgroups={len(g.lattice())} sd=sdStar={g.sd()} '
              'center=1 proper_subgroups_abelian=True normals=3')

    # Coprime subgroup splitting and section-minimum product, direct small check.
    assert a4 is not None
    c5 = Group(range(5), lambda x, y: (x + y) % 5)
    g = Group(product(range(a4.n), range(5)),
              lambda x, y: (a4.t[x[0]][y[0]], (x[1] + y[1]) % 5))
    assert gcd(a4.n, c5.n) == 1
    assert len(g.lattice()) == len(a4.lattice()) * len(c5.lattice())
    for h in g.lattice():
        coords = [g.elements[x] for x in h]
        assert set(coords) == set(product({x for x, _ in coords}, {y for _, y in coords}))
    assert g.sd() == a4.sd() * c5.sd()
    assert min(g.section_values()) == min(a4.section_values()) * min(c5.section_values())
    print('ENUMERATED_COPRIME A4 x C5: subgroups=20 sd=sdStar=16/25')

    # Explicit large witness: exact arithmetic on the classified subgroup lattice.
    assert prime(127) and order(2, 127) == 7
    assert gaussian_row(2, 7) == [1, 127, 2667, 11811, 11811, 2667, 127, 1]
    assert sum(gaussian_row(2, 7)) == 29212
    s = scalar_sd(2, 7)
    assert s == F(853400265, 860894281) and s > F(99, 100)
    assert 100 * s.numerator - 99 * s.denominator == 111492681
    assert [pow(31, i, 11) for i in range(1, 6)] == [9, 4, 3, 5, 1]
    assert order(31, 11) == 5
    a = sum(gaussian_row(31, 5))
    assert a == 1837991432
    assert a == 2*31**6 + 2*31**5 + 6*31**4 + 6*31**3 + 6*31**2 + 4*31 + 6
    t = scalar_sd(31, 5)
    assert t == F(1689106253990898773, 1742136202306250528)
    assert t > F(62, 63)**2 > F(24, 25)
    assert 25 * 62**2 - 24 * 63**2 == 844
    assert gcd(2**7 * 127, 31**5 * 11) == 1
    assert s*t > F(2376, 2500) > F(23, 25)
    print(f'EXPLICIT S: a=29212 m=128 sdStar_formula={s}')
    print(f'EXPLICIT T: a={a} m={31**5} sdStar_formula={t}')
    print(f'EXPLICIT product: {s*t} > 2376/2500 > 23/25')

    # Finite family samples. Uniform existence is supplied by the prose argument,
    # not claimed as a consequence of running a finite program.
    for c, k in [(F(23,25),2), (F(99,100),3), (F(9999,10000),5)]:
        ps = [2,3,5,7,11][:k]
        used = set(ps)
        lower_product, parameters = F(1), []
        target = 1 - (1-c)/k
        for p in ps:
            for q in range(2, 2000):
                if q in used or not prime(q):
                    continue
                r = order(p,q)
                if r <= 4:
                    continue
                bound_a = p ** (r*r//4)
                lower = F(bound_a+1, bound_a+1+p**r)**2
                if lower > target:
                    used.add(q)
                    lower_product *= lower
                    parameters.append((p,q,r))
                    break
            else:
                raise AssertionError('finite sample search exhausted')
        assert lower_product > c
        assert all(gcd(p**r*q, u**v*w) == 1 for i,(p,q,r) in enumerate(parameters)
                   for u,w,v in parameters[i+1:])
        print(f'FAMILY_SAMPLE c={c} k={k} (p,q,r)={parameters} certified_product_lower_gt_c=True')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
