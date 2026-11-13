#!/usr/bin/env python3
"""Exact Haar MST checks. Only integer and Fraction arithmetic; no simulation."""
from fractions import Fraction as Q
from itertools import combinations, product
from math import factorial, isqrt


def sqrt_split(n):
    a, b, p = 1, n, 2
    while p * p <= b:
        while b % (p * p) == 0:
            a *= p
            b //= p * p
        p += 1
    return a, b


class Radical:
    """Exact rational linear combinations of square roots of square-free integers."""
    def __init__(self, x=0):
        c = dict(x) if isinstance(x, dict) else ({1: Q(x)} if x else {})
        self.c = {b: Q(v) for b, v in c.items() if v}

    @staticmethod
    def lift(x):
        return x if isinstance(x, Radical) else Radical(x)

    def __add__(self, x):
        out = dict(self.c)
        for b, c in self.lift(x).c.items():
            out[b] = out.get(b, Q(0)) + c
        return Radical(out)
    __radd__ = __add__

    def __neg__(self):
        return Radical({b: -c for b, c in self.c.items()})

    def __sub__(self, x):
        return self + -self.lift(x)

    def __mul__(self, x):
        out = {}
        for b, c in self.c.items():
            for d, e in self.lift(x).c.items():
                a, f = sqrt_split(b * d)
                out[f] = out.get(f, Q(0)) + a * c * e
        return Radical(out)
    __rmul__ = __mul__

    def __truediv__(self, x):
        return self * (Q(1) / x)

    def __eq__(self, x):
        return self.c == self.lift(x).c


def frequency(k):
    a, b = sqrt_split(k)
    return Radical({b: k * a})


def source_checks():
    for d in range(2, 13):
        lam = [frequency(k) for k in range(1, d + 1)]
        squares = sum((x * x for x in lam), Radical())
        total = sum(lam, Radical())
        variance = squares / d - total * total / (d * d)
        # Independently sum Dirichlet(1,...,1) second moments of the weighted mean.
        mean_square = sum((lam[i] * lam[j] * Q(1 + (i == j), d * (d + 1))
                           for i in range(d) for j in range(d)), Radical())
        marginal_second = (squares / d - mean_square) / (d - 1)
        assert marginal_second == variance * Q(d, (d - 1) * (d + 1))
        # Actual rational orthogonal conjugations, with the source's exact k^(3/2).
        norm2 = sum(k * k for k in range(1, d + 1))
        u = [[Q(i == j) - Q(2 * (i + 1) * (j + 1), norm2)
              for j in range(d)] for i in range(d)]
        for j in range(d):
            h = [sum((lam[a] * u[a][i] * u[a][j] for a in range(d)), Radical())
                 for i in range(d)]
            lhs = sum((h[i] * h[i] for i in range(d) if i != j), Radical())
            m1 = sum((lam[a] * u[a][j] * u[a][j] for a in range(d)), Radical())
            m2 = sum((lam[a] * lam[a] * u[a][j] * u[a][j]
                      for a in range(d)), Radical())
            assert lhs == m2 - m1 * m1
    # Integral_0^1 x^p dx = 1/(p+1), p>-1.
    m1, m2 = Q(1) / (Q(3, 2) + 1), Q(1) / (Q(3) + 1)
    v = m2 - m1 * m1
    assert v == Q(9, 100)
    r = Q(isqrt(v.numerator), isqrt(v.denominator))
    assert r * r == v and 1 / r == Q(10, 3)
    print('FIDELITY: exact source marginal second moment for every d=2..12')
    print('FIDELITY: spectral variance coefficient 1/4-(2/5)^2=9/100')
    print('CONSTANT: reciprocal square root of 9/100 equals 10/3 exactly')
    print('COLUMN ALGEBRA: exact rational orthogonal conjugations, every column, d=2..12')


def radius_checks():
    # Dirichlet grouped inverse moment, only when group size m>=2.
    for d in range(16, 1025):
        m = d // 4
        inverse = Q((d - 1) * (d - 2), (m - 1) ** 2)
        assert m >= 2 and 4 * inverse <= 256
    # (3 sqrt(3)-1)/8 > 1/2 follows from 27>25.
    assert 27 > 25
    print('INVERSE RADIUS: grouped Dirichlet bound <=256 at every d=16..1024')


def prufer_tree(n, code):
    degrees = [1] * n
    for x in code:
        degrees[x] += 1
    edges = []
    for x in code:
        leaf = next(i for i in range(n) if degrees[i] == 1)
        edges.append(tuple(sorted((leaf, x))))
        degrees[leaf] -= 1
        degrees[x] -= 1
    a, b = [i for i in range(n) if degrees[i] == 1]
    return edges + [(a, b)]


def is_tree(n, edges):
    if len(edges) != n - 1:
        return False
    parent = list(range(n))
    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i
    for a, b in edges:
        ra, rb = root(a), root(b)
        if ra == rb:
            return False
        parent[ra] = rb
    return len({root(i) for i in range(n)}) == 1


def mst_checks():
    checked = 0
    for n in range(2, 6):
        edges = list(combinations(range(n), 2))
        trees = [prufer_tree(n, p) for p in product(range(n), repeat=n - 2)]
        assert len(trees) == n ** (n - 2) and all(is_tree(n, t) for t in trees)
        alphabet = (1, 2, 3) if n <= 4 else (1, 2)
        for costs in product(alphabet, repeat=len(edges)):
            w = dict(zip(edges, costs))
            minimum = min(sum(w[e] for e in t) for t in trees)
            for k in range(1, n):
                t = [(i, i + 1) for i in range(k - 1)]
                for j in range(k, n):
                    a = min(range(k), key=lambda a: (w[(a, j)], a))
                    t.append((a, j))
                assert is_tree(n, t)
                assert minimum <= sum(w[e] for e in t)
            for threshold in (Q(1, 2), Q(3, 2), Q(5, 2), Q(7, 2)):
                cheap = sum(c < threshold for c in costs)
                assert minimum >= threshold * (n - 1 - cheap)
            checked += 1
    print(f'ANCHOR UPPER / CHEAP-EDGE LOWER: {checked} weighted complete graphs exhaustively checked')
    print('DOMAINS: n=2..4 weights={1,2,3}; n=5 weights={1,2}; every anchor size; four thresholds')


def main():
    source_checks()
    radius_checks()
    mst_checks()
    print('LIMIT: exact finite identities and tree inequalities; no asymptotic proof or Haar simulation')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
