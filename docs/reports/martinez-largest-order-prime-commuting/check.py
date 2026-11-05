#!/usr/bin/env python3
"""Exact finite-group checks for Martínez's largest-order-prime question.

Standard-library only. Group tables, orders, normal closures, and fractions are
computed from multiplication, without substituting the proposed count formulae.
The finite checks support a separate general argument; they do not prove it.
"""
from collections import Counter, deque
from fractions import Fraction
from itertools import permutations, product


def prime_support(n):
    assert n >= 1
    result = set()
    d = 2
    while d * d <= n:
        if n % d == 0:
            result.add(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        result.add(n)
    return result


class Group:
    def __init__(self, elements, multiply, identity, check_associativity=True):
        self.elements = tuple(elements)
        self.index = {x: i for i, x in enumerate(self.elements)}
        assert len(self.index) == len(self.elements)
        self.e = self.index[identity]
        self.n = len(self.elements)
        self.table = [[self.index[multiply(x, y)] for y in self.elements]
                      for x in self.elements]
        t, e, n = self.table, self.e, self.n
        assert all(t[e][a] == a == t[a][e] for a in range(n))
        self.inv = [next(b for b in range(n) if t[a][b] == e == t[b][a])
                    for a in range(n)]
        if check_associativity:
            assert all(t[t[a][b]][c] == t[a][t[b][c]]
                       for a in range(n) for b in range(n) for c in range(n))
        self.orders = []
        for a in range(n):
            v = e
            for k in range(1, n + 1):
                v = t[v][a]
                if v == e:
                    self.orders.append(k)
                    break
            else:
                raise AssertionError('element order exceeds group size')
        self.normal_closures = {}

    def pi_elements(self, pi):
        return {a for a, order in enumerate(self.orders)
                if prime_support(order) <= pi}

    def centralizer(self, x):
        t = self.table
        return {a for a in range(self.n) if t[a][x] == t[x][a]}

    def center(self):
        return {a for a in range(self.n) if len(self.centralizer(a)) == self.n}

    def normal_closure(self, x):
        t = self.table
        generators = frozenset(t[t[g][x]][self.inv[g]] for g in range(self.n))
        if generators not in self.normal_closures:
            seen = {self.e}
            queue = deque([self.e])
            while queue:
                a = queue.popleft()
                for b in generators:
                    c = t[a][b]
                    if c not in seen:
                        seen.add(c)
                        queue.append(c)
            self.normal_closures[generators] = frozenset(seen)
        return self.normal_closures[generators]

    def pi_core(self, pi):
        # x is in a normal pi-subgroup iff its normal closure is a pi-group.
        core = {a for a in range(self.n)
                if prime_support(len(self.normal_closure(a))) <= pi}
        t = self.table
        assert self.e in core
        assert all(self.inv[a] in core for a in core)
        assert all(t[a][b] in core for a in core for b in core)
        assert all(t[t[g][a]][self.inv[g]] in core
                   for a in core for g in range(self.n))
        assert prime_support(len(core)) <= pi
        return core

    def commuting_fraction(self, pi, x):
        pe = self.pi_elements(pi)
        assert self.e in pe
        return Fraction(len(self.centralizer(x) & pe), len(pe))


def symmetric_three():
    return Group(permutations(range(3)),
                 lambda a, b: tuple(a[b[i]] for i in range(3)), (0, 1, 2))


def dihedral(q):
    return Group(product(range(q), range(2)),
                 lambda a, b: ((a[0] + (-1 if a[1] else 1) * b[0]) % q,
                               (a[1] + b[1]) % 2), (0, 0))


def direct_product(h, k):
    # A direct product of the two exhaustively checked associative factors.
    return Group(product(h.elements, k.elements),
                 lambda a, b: (h.elements[h.table[h.index[a[0]]][h.index[b[0]]]],
                               k.elements[k.table[k.index[a[1]]][k.index[b[1]]]]),
                 (h.elements[h.e], k.elements[k.e]), check_associativity=False)


def psl2_five():
    def canonical(a):
        return min(a, tuple((-x) % 5 for x in a))

    elements = sorted({canonical((a, b, c, d))
                       for a, b, c, d in product(range(5), repeat=4)
                       if (a * d - b * c) % 5 == 1})

    def multiply(a, b):
        return canonical(((a[0] * b[0] + a[1] * b[2]) % 5,
                          (a[0] * b[1] + a[1] * b[3]) % 5,
                          (a[2] * b[0] + a[3] * b[2]) % 5,
                          (a[2] * b[1] + a[3] * b[3]) % 5))
    return Group(elements, multiply, (1, 0, 0, 1))


def main():
    # Source, introduction p. 3: Pr_p(PSL(2,p))=(p^2+p-1)/p^3 for p>3.
    g = psl2_five()
    pe = g.pi_elements({5})
    pairs = sum(len(g.centralizer(x) & pe) for x in pe)
    probability = Fraction(pairs, len(pe) ** 2)
    assert g.n == 60 and len(pe) == 25 and pairs == 145
    assert probability == Fraction(5 ** 2 + 5 - 1, 5 ** 3) == Fraction(29, 125)
    print(f'FIDELITY PSL(2,5): order={g.n} pi_elements={len(pe)} '
          f'commuting_pairs={pairs} Pr_5={probability}')

    s3 = symmetric_three()
    tau = (1, 0, 2)
    t = s3.index[tau]
    assert len(s3.pi_elements({2})) == 4
    assert s3.pi_core({2}) == {s3.e}
    assert s3.commuting_fraction({2}, t) == Fraction(1, 2)
    assert s3.center() == {s3.e}
    assert s3.pi_elements(set()) == {s3.e}
    print('FIDELITY Theorem A, S3 and pi={2}: fraction=1/2; pi_core=1')

    for q in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        assert prime_support(q) == {q}
        d = dihedral(q)
        r = d.index[(1, 0)]
        pi = {2, q}
        assert d.pi_elements(pi) == set(range(d.n))
        assert len(d.centralizer(r)) == q and d.center() == {d.e}
        g = direct_product(s3, d)
        x = g.index[(tau, (1, 0))]
        pe, core = g.pi_elements(pi), g.pi_core(pi)
        central = g.centralizer(x) & pe
        expected_core = {g.index[(s3.elements[s3.e], y)] for y in d.elements}
        assert core == expected_core and x not in core and x in pe
        assert g.n == 12 * q and g.orders[x] == 2 * q
        assert max(prime_support(g.orders[x])) == q
        assert len(pe) == 8 * q and len(central) == 2 * q
        assert g.center() == {g.e}
        assert len(g.normal_closure(x)) == 6 * q
        ratio = g.commuting_fraction(pi, x)
        assert ratio == Fraction(1, 4) > Fraction(1, q)
        assert all(g.commuting_fraction(pi, y) <= Fraction(1, 2)
                   for y in pe - core)
        print(f'COUNTEREXAMPLE q={q}: order={g.n} order_x={g.orders[x]} '
              f'pi_elements={len(pe)} centralizing_pi_elements={len(central)} '
              f'pi_core={len(core)} normal_closure_x={len(g.normal_closure(x))} '
              f'center=1 ratio={ratio} > 1/{q}; Theorem_A_OK')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
