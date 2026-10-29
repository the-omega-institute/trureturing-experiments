#!/usr/bin/env python3
"""Exact finite checks for ALM arXiv:2209.00986v2, Problem 5.2.

No external dependencies. Probabilities are computed by a weighted permanent
of the actual left/right coset intersection matrix, independently of the
published factorial product. Every subgroup and every normal subgroup of
each explicitly constructed group is enumerated.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
from math import factorial, prod


class Group:
    def __init__(self, name, elements, mul, identity):
        self.name = name
        self.elements = list(elements)
        self.n = len(self.elements)
        position = {x: i for i, x in enumerate(self.elements)}
        self.one = position[identity]
        self.table = [[position[mul(a, b)] for b in self.elements]
                      for a in self.elements]
        self.inv = [next(b for b in range(self.n)
                         if self.table[a][b] == self.one
                         and self.table[b][a] == self.one)
                    for a in range(self.n)]
        m = self.table
        assert all(m[self.one][a] == a == m[a][self.one]
                   for a in range(self.n))
        assert all(m[m[a][b]][c] == m[a][m[b][c]]
                   for a, b, c in product(range(self.n), repeat=3))

    def closure(self, seeds):
        generators = tuple(seeds)
        found, queue = {self.one}, [self.one]
        for a in queue:
            for b in generators:
                c = self.table[a][b]
                if c not in found:
                    found.add(c)
                    queue.append(c)
        return frozenset(found)

    @lru_cache(None)
    def subgroups(self):
        seen = {frozenset([self.one])}
        queue = list(seen)
        for h in queue:
            for g in range(self.n):
                if g not in h:
                    k = self.closure(tuple(sorted(h)) + (g,))
                    if k not in seen:
                        seen.add(k)
                        queue.append(k)
        return tuple(sorted(seen, key=lambda h: (len(h), tuple(sorted(h)))))

    def normal(self, h):
        return all(self.table[self.table[g][a]][self.inv[g]] in h
                   for g in range(self.n) for a in h)

    @lru_cache(None)
    def cosets(self, h, left=True):
        unseen, result = set(range(self.n)), []
        while unseen:
            g = min(unseen)
            c = frozenset(self.table[g][a] if left else self.table[a][g]
                          for a in h)
            assert c <= unseen
            unseen -= c
            result.append(c)
        return tuple(result)

    @lru_cache(None)
    def probability(self, h):
        left, right = self.cosets(h), self.cosets(h, False)
        rows = tuple(tuple((j, len(c & d)) for j, d in enumerate(right) if c & d)
                     for c in left)
        # Each chosen cell contributes one element. Distinct right-coset
        # columns enforce that each right coset is met exactly once.
        @lru_cache(None)
        def permanent(i, used):
            if i == len(rows):
                return 1
            return sum(weight * permanent(i + 1, used | (1 << j))
                       for j, weight in rows[i] if not used & (1 << j))
        return Fraction(permanent(0, 0), len(h) ** len(left))

    def orbits(self, acting, stabilizer):
        cosets = self.cosets(stabilizer)
        index = {g: i for i, c in enumerate(cosets) for g in c}
        unseen, result = set(range(len(cosets))), []
        while unseen:
            i = min(unseen)
            g = min(cosets[i])
            orbit = frozenset(index[self.table[h][g]] for h in acting)
            assert orbit <= unseen
            unseen -= orbit
            result.append(orbit)
        return tuple(result)

    def quotient(self, normal):
        assert self.normal(normal)
        cosets = self.cosets(normal)
        index = {g: i for i, c in enumerate(cosets) for g in c}
        representatives = [min(c) for c in cosets]
        q = Group(self.name + '/N', range(len(cosets)),
                  lambda i, j: index[self.table[representatives[i]][representatives[j]]],
                  index[self.one])
        return q, index


def cyclic(n):
    return Group('C' + str(n), range(n), lambda a, b: (a + b) % n, 0)


def abelian(*orders):
    return Group('x'.join('C' + str(n) for n in orders),
                 product(*(range(n) for n in orders)),
                 lambda a, b: tuple((x + y) % n for x, y, n in zip(a, b, orders)),
                 (0,) * len(orders))


def dihedral(n):
    return Group('D' + str(n) + '(order ' + str(2*n) + ')',
                 product(range(n), range(2)),
                 lambda a, b: ((a[0] + (-1) ** a[1] * b[0]) % n,
                               (a[1] + b[1]) % 2), (0, 0))


def quaternion():
    def mul(a, b):
        s, i = a
        t, j = b
        if i == 0:
            return s*t, j
        if j == 0:
            return s*t, i
        if i == j:
            return -s*t, 0
        return s*t*(1 if (i, j) in [(1, 2), (2, 3), (3, 1)] else -1), 6-i-j
    return Group('Q8', product((1, -1), range(4)), mul, (1, 0))


def permutation_group(n, even=False):
    elements = list(permutations(range(n)))
    if even:
        elements = [p for p in elements if
                    sum(p[i] > p[j] for i in range(n) for j in range(i+1, n)) % 2 == 0]
    return Group(('A' if even else 'S') + str(n), elements,
                 lambda a, b: tuple(a[b[i]] for i in range(n)), tuple(range(n)))


def direct_product(a, b):
    return Group(a.name + 'x' + b.name, product(range(a.n), range(b.n)),
                 lambda x, y: (a.table[x[0]][y[0]], b.table[x[1]][y[1]]),
                 (a.one, b.one))


def dicyclic12():
    # a^6=1, b^2=a^3, b*a*b^-1=a^-1; element a^i*b^j.
    return Group('Dic12', product(range(6), range(2)),
                 lambda a, b: ((a[0] + (-1)**a[1]*b[0] + 3*a[1]*b[1]) % 6,
                               (a[1] + b[1]) % 2), (0, 0))


def a(t):
    assert t >= 1
    return Fraction(factorial(t), t**t)


def check_group(g):
    subgroups = g.subgroups()
    normals = [n for n in subgroups if g.normal(n)]
    subgroup_set = set(subgroups)
    for h in subgroups:
        # Source Proposition 3.2, checked against the independent permanent.
        assert g.probability(h) == prod((a(len(o)) for o in g.orbits(h, h)), start=Fraction(1))
        assert (g.probability(h) == 1) == g.normal(h)  # Source Corollary 3.4.
    comparisons = strict = 0
    for n in normals:
        quotient, quotient_map = g.quotient(n)
        for h in subgroups:
            l = frozenset(g.table[x][y] for x in h for y in n)
            assert l in subgroup_set
            ph, pl = g.probability(h), g.probability(l)
            assert ph <= pl
            assert (ph == pl) == (n <= h or g.normal(h))
            image_h = frozenset(quotient_map[x] for x in h)
            assert quotient.probability(image_h) == pl
            # Verify the specific orbit/fiber bridge on actual cosets.
            upstairs, downstairs = g.cosets(h), g.cosets(l)
            down_index = {x: i for i, c in enumerate(downstairs) for x in c}
            pi = [down_index[min(c)] for c in upstairs]
            q = len(l) // len(h)
            assert all(len({down_index[x] for x in c}) == 1 for c in upstairs)
            assert all(pi.count(j) == q for j in range(len(downstairs)))
            down_orbits = g.orbits(h, l)
            assert set(down_orbits) == set(g.orbits(l, l))
            up_orbits = g.orbits(h, h)
            for b in down_orbits:
                pieces = [o for o in up_orbits if pi[min(o)] in b]
                assert pieces
                assert all(frozenset(pi[i] for i in o) == b for o in pieces)
                assert sum(map(len, pieces)) == q * len(b)
                lhs = prod((a(len(o)) for o in pieces), start=Fraction(1))
                assert lhs <= a(len(b))
                if q > 1 and len(b) > 1:
                    assert lhs < a(len(b))
            comparisons += 1
            strict += ph < pl
    print(f'{g.name}: order={g.n} subgroups={len(subgroups)} normals={len(normals)} '
          f'comparisons={comparisons} strict={strict} tp={min(map(g.probability, subgroups))}')
    return comparisons, strict


def main():
    # The standard list of all 24 isomorphism types of order at most 12.
    groups = [cyclic(n) for n in range(1, 13)]
    groups += [abelian(2, 2), dihedral(3), abelian(4, 2), abelian(2, 2, 2),
               dihedral(4), quaternion(), abelian(3, 3), dihedral(5),
               abelian(6, 2), dihedral(6), permutation_group(4, even=True), dicyclic12()]
    assert len(groups) == 24
    # Additional nonabelian groups give nontrivial quotient orbit patterns.
    groups += [dihedral(7), dihedral(8), permutation_group(4),
               direct_product(permutation_group(3), permutation_group(3))]
    total = strict = 0
    for g in groups:
        c, s = check_group(g)
        total += c
        strict += s
    # Fidelity: source Lemma 4.2(ii), D_n has order 2n.
    for n in range(3, 9):
        g = next(g for g in groups if g.name == f'D{n}(order {2*n})')
        expected = Fraction(1, 2) ** ((n-1)//2 if n % 2 else (n-2)//2)
        actual = min(map(g.probability, g.subgroups()))
        assert actual == expected
        print(f'FIDELITY Lemma 4.2(ii): tp(D_{n})={actual}')
    g = next(g for g in groups if g.name == 'A4')
    assert min(map(g.probability, g.subgroups())) == Fraction(2, 9)
    print('FIDELITY Lemma 4.2(iv): tp(A4)=2/9')
    for t in range(1, 101):
        assert a(t+1)/a(t) == Fraction(t, t+1)**t < 1
    print(f'SUMMARY groups={len(groups)} all_order_at_most_12_types=24 '
          f'subgroup_normal_pairs={total} strict_pairs={strict}')
    print('All arithmetic evidence uses integers and fractions.Fraction; no floating point.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
