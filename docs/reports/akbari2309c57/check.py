"""Exact small finite-group checks for arXiv:2309.09104v3, Conjecture 5.7."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import lcm


class GroupTable:
    """Multiplication, generation, and derived-series operations on a finite group."""

    def __init__(self, elements, multiply, identity):
        self.elements = tuple(elements)
        self.index = {element: i for i, element in enumerate(self.elements)}
        self.one = self.index[identity]
        self.indices = range(len(self.elements))
        self.mul = [
            [self.index[multiply(x, y)] for y in self.elements]
            for x in self.elements
        ]
        self.inverse = [row.index(self.one) for row in self.mul]
        self.solvable = {}

    def generated(self, generators):
        generators = tuple(set(generators) - {self.one})
        seen = {self.one}
        pending = [self.one]
        while pending:
            row = self.mul[pending.pop()]
            for g in generators:
                h = row[g]
                if h not in seen:
                    seen.add(h)
                    pending.append(h)
        return frozenset(seen)

    def derived(self, subgroup):
        commutators = {
            self.mul[self.mul[self.mul[x][y]][self.inverse[x]]][self.inverse[y]]
            for x, y in combinations(subgroup, 2)
        }
        return self.generated(commutators)

    def is_soluble(self, subgroup):
        if subgroup not in self.solvable:
            if len(subgroup) == 1:
                self.solvable[subgroup] = True
            else:
                derived = self.derived(subgroup)
                self.solvable[subgroup] = (
                    derived != subgroup and self.is_soluble(derived)
                )
        return self.solvable[subgroup]

    def solubilizer(self, x):
        return {
            y for y in self.indices if self.is_soluble(self.generated((x, y)))
        }

    def radical(self):
        # An element lies in the radical exactly when its normal closure is soluble.
        return {
            x
            for x in self.indices
            if self.is_soluble(
                self.generated(
                    self.mul[self.mul[g][x]][self.inverse[g]] for g in self.indices
                )
            )
        }

    def order(self, x):
        return len(self.generated((x,)))

    def is_insoluble(self):
        return not self.is_soluble(frozenset(self.indices))


def perm_mul(x, y):
    return tuple(x[y_i] for y_i in y)


def even_permutations(n):
    return [
        p for p in permutations(range(n))
        if sum(p[i] > p[j] for i, j in combinations(range(n), 2)) % 2 == 0
    ]


def check_a5():
    identity = tuple(range(5))
    a = (1, 0, 3, 2, 4)
    group = GroupTable(even_permutations(5), perm_mul, identity)
    assert len(group.elements) == 60
    assert group.is_insoluble()
    assert group.radical() == {group.one}
    counts = Counter()
    for x in group.indices:
        size = len(group.solubilizer(x))
        order = group.order(x)
        fraction = Fraction(size, 60)
        assert (fraction > Fraction(1, 2)) == (order in (1, 2))
        if order == 2:
            assert fraction == Fraction(3, 5)
        counts[order, size] += 1
    assert counts == {(1, 60): 1, (2, 36): 15, (3, 24): 20, (5, 10): 24}
    print("Corollary 5.6 A5 row, (order, solubilizer size) -> multiplicity:")
    print(sorted(counts.items()))

    soluble_partners = {group.elements[y] for y in group.solubilizer(group.index[a])}
    by_order = Counter(group.order(group.index[y]) for y in soluble_partners)
    assert by_order == {1: 1, 2: 15, 3: 12, 5: 8}
    print("A5 witness, soluble partners by order:", sorted(by_order.items()))
    return group.elements, identity, a, soluble_partners


def check_products(perms, identity, a, soluble_partners):
    for m in (3, 5, 7, 9):
        group = GroupTable(
            product(perms, range(m)),
            lambda x, y: (perm_mul(x[0], y[0]), (x[1] + y[1]) % m),
            (identity, 0),
        )
        x = group.index[a, 1]
        radical = group.radical()
        solubilizer = group.solubilizer(x)
        assert len(group.elements) == 60 * m and group.is_insoluble()
        assert group.order(x) == lcm(2, m) == 2 * m
        assert x not in radical
        assert {group.elements[y] for y in radical} == set(product((identity,), range(m)))
        assert {group.elements[y] for y in solubilizer} == set(
            product(soluble_partners, range(m))
        )
        assert len(solubilizer) == 36 * m
        probability = Fraction(len(solubilizer), len(group.elements))
        assert probability == Fraction(3, 5) > Fraction(1, 2)
        print(
            f"A5 x C{m}: |G|={len(group.elements)} order(x)={group.order(x)} "
            f"|R(G)|={len(radical)} |Sol_G(x)|={len(solubilizer)} P_S(x)={probability}"
        )


def matrix_mul(x, y):
    a, b, c, d = x
    e, f, g, h = y
    return ((a * e + b * g) % 5, (a * f + b * h) % 5,
            (c * e + d * g) % 5, (c * f + d * h) % 5)


def check_sl2():
    matrices = [
        (a, b, c, d) for a, b, c, d in product(range(5), repeat=4)
        if (a * d - b * c) % 5 == 1
    ]
    group = GroupTable(matrices, matrix_mul, (1, 0, 0, 1))
    x = group.index[0, 4, 1, 0]
    radical = group.radical()
    solubilizer = group.solubilizer(x)
    assert len(group.elements) == 120 and group.is_insoluble()
    assert group.order(x) == 4 and x not in radical
    assert group.mul[x][x] == group.index[4, 0, 0, 4]
    assert {group.elements[y] for y in radical} == {(1, 0, 0, 1), (4, 0, 0, 4)}
    assert len(solubilizer) == 72
    probability = Fraction(len(solubilizer), len(group.elements))
    assert probability == Fraction(3, 5) > Fraction(1, 2)
    print(
        f"SL2(F5): |G|={len(group.elements)} order(x)={group.order(x)} "
        f"|R(G)|={len(radical)} |Sol_G(x)|={len(solubilizer)} P_S(x)={probability}"
    )


if __name__ == "__main__":
    check_products(*check_a5())
    check_sl2()
    print("ALL_CHECKS_OK")
