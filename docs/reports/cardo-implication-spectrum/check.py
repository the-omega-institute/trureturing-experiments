#!/usr/bin/env python3
"""Exact finite checks for Cardo's native implication-spectrum question.

Python integers encode truth tables; bit a is the value at valuation a.
No external packages, random choices, or floating-point arithmetic are used.
These checks do not prove the unbounded spectrum assertion.
"""

from fractions import Fraction


def coordinate(k, i):
    return sum(1 << a for a in range(1 << k) if (a >> i) & 1)


def closure(k):
    """Least truth-table set containing projections and closed under implication."""
    full = (1 << (1 << k)) - 1
    tables = set(coordinate(k, i) for i in range(k))
    pending = list(tables)
    while pending:
        a = pending.pop()
        for b in tuple(tables):
            for c in ((full ^ a) | b, (full ^ b) | a):
                if c not in tables:
                    tables.add(c)
                    pending.append(c)
    return tables


class NativeTerms:
    """Interned acyclic syntax with exactly variable and binary-implication nodes."""

    def __init__(self, k):
        self.k = k
        self.full = (1 << (1 << k)) - 1
        self.nodes = []
        self.tables = []
        self.supports = []
        self.intern = {}
        self.variables = [self.var(i) for i in range(k)]

    def var(self, i):
        node = ("var", i)
        if node not in self.intern:
            self.intern[node] = len(self.nodes)
            self.nodes.append(node)
            self.tables.append(coordinate(self.k, i))
            self.supports.append(1 << i)
        return self.intern[node]

    def imp(self, a, b):
        node = ("imp", a, b)
        if node not in self.intern:
            self.intern[node] = len(self.nodes)
            self.nodes.append(node)
            self.tables.append((self.full ^ self.tables[a]) | self.tables[b])
            self.supports.append(self.supports[a] | self.supports[b])
        return self.intern[node]

    def join(self, a, b):
        return self.imp(self.imp(a, b), b)

    def guarded_synthesis(self, target, guard):
        """DNF with temporary zero replaced by z, followed by z OR phi."""
        z = self.variables[guard]
        assert (target & coordinate(self.k, guard)) == coordinate(self.k, guard)
        one = self.imp(z, z)
        phi = z
        for valuation in range(1 << self.k):
            if ((valuation >> guard) & 1) or not ((target >> valuation) & 1):
                continue
            conjunction = one
            for i in range(self.k):
                if i == guard:
                    continue
                x = self.variables[i]
                literal = x if (valuation >> i) & 1 else self.imp(x, z)
                conjunction = self.imp(self.imp(conjunction, self.imp(literal, z)), z)
            phi = self.join(phi, conjunction)
        return self.imp(self.imp(z, phi), phi)

    def use_all_variables(self, t):
        for i, x in enumerate(self.variables):
            if not ((self.supports[t] >> i) & 1):
                t = self.imp(self.imp(x, x), t)
        assert self.supports[t] == (1 << self.k) - 1
        return t

    def scalar_cross_check(self):
        """Re-evaluate every syntax node at every valuation using scalar booleans."""
        for valuation in range(1 << self.k):
            values = []
            for j, node in enumerate(self.nodes):
                if node[0] == "var":
                    value = bool((valuation >> node[1]) & 1)
                else:
                    assert node[0] == "imp"
                    a, b = node[1:]
                    assert a < j and b < j
                    value = (not values[a]) or values[b]
                values.append(value)
                assert value == bool((self.tables[j] >> valuation) & 1)


def fidelity():
    """Source Section 4's semilattice formula, with nonempty y and z blocks."""
    cases = 0
    for p in range(1, 4):
        for q in range(1, 4):
            for r in range(3):
                k = p + q + r
                count = 0
                for a in range(1 << k):
                    x = all((a >> i) & 1 for i in range(r))
                    y = all((a >> i) & 1 for i in range(r, r + p))
                    z = all((a >> i) & 1 for i in range(r + p, k))
                    count += ((x and y) == (x and z))
                observed = Fraction(count, 1 << k)
                expected = 1 - Fraction((1 << p) + (1 << q) - 2, 1 << k)
                assert observed == expected
                cases += 1
    assert Fraction(6, 8) == Fraction(3, 4)
    print(f"FIDELITY semilattice phi(p,q,r): {cases} exact cases; p,q=1..3, r=0..2")


def all_small_term_functions():
    for k in range(1, 5):
        n = 1 << k
        full = (1 << n) - 1
        tables = closure(k)
        guarded = set()
        syntax = NativeTerms(k)
        syntheses = 0
        for guard in range(k):
            low = [a for a in range(n) if not ((a >> guard) & 1)]
            for subset in range(1 << len(low)):
                zeros = sum(1 << a for i, a in enumerate(low) if (subset >> i) & 1)
                target = full ^ zeros
                guarded.add(target)
                term = syntax.use_all_variables(syntax.guarded_synthesis(target, guard))
                assert syntax.tables[term] == target
                syntheses += 1
        assert tables == guarded
        counts = {n - (s ^ t).bit_count() for s in tables for t in tables}
        expected_counts = set(range((n + 3) // 4, n + 1))
        assert counts == expected_counts
        syntax.scalar_cross_check()
        probabilities = {Fraction(c, n) for c in counts}
        print(f"CLOSURE k={k}: {len(tables)} term functions; {syntheses} guarded syntheses; "
              f"{len(counts)} agreement counts; min={min(probabilities)}, max={max(probabilities)}")


def every_small_dyadic():
    total = 0
    for k in range(2, 7):
        n = 1 << k
        full = (1 << n) - 1
        union = [a for a in range(n) if not (a & 1) or not ((a >> 1) & 1)]
        assert len(union) == 3 * (1 << (k - 2))
        syntax = NativeTerms(k)
        cases = 0
        for d in range(len(union) + 1):
            disagreement = union[:d]
            d1 = [a for a in disagreement if not (a & 1)]
            d2 = [a for a in disagreement if a & 1]
            assert all(not ((a >> 1) & 1) for a in d2)
            assert set(d1).isdisjoint(d2) and set(d1) | set(d2) == set(disagreement)
            f1 = full ^ sum(1 << a for a in d1)
            f2 = full ^ sum(1 << a for a in d2)
            s = syntax.use_all_variables(syntax.guarded_synthesis(f1, 0))
            t = syntax.use_all_variables(syntax.guarded_synthesis(f2, 1))
            assert syntax.tables[s] == f1 and syntax.tables[t] == f2
            observed = Fraction(n - (syntax.tables[s] ^ syntax.tables[t]).bit_count(), n)
            assert observed == Fraction(n - d, n)
            assert Fraction(1, 4) <= observed <= 1
            cases += 1
        syntax.scalar_cross_check()
        total += cases
        print(f"REALIZATION k={k}: every numerator {n // 4}..{n} over denominator {n}; "
              f"{cases} native equation pairs using all {k} variables")
    print(f"REALIZATION_TOTAL {total}")


def endpoints():
    syntax = NativeTerms(2)
    x, y = syntax.variables
    t = syntax.imp(x, y)
    assert Fraction(4 - (syntax.tables[x] ^ syntax.tables[t]).bit_count(), 4) == Fraction(1, 4)
    assert Fraction(4 - (syntax.tables[x] ^ syntax.tables[x]).bit_count(), 4) == 1
    syntax.scalar_cross_check()
    print("ENDPOINTS x = (x -> y): 1/4; x = x: 1")


if __name__ == "__main__":
    fidelity()
    all_small_term_functions()
    every_small_dyadic()
    endpoints()
    print("ALL_CHECKS_OK")
