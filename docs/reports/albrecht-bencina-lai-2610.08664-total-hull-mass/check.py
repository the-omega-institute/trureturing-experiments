#!/usr/bin/env python3
"""Exact finite-field checks for the proposed total-Galois-hull mass formula.

Standard-library only. Subspaces are enumerated once via their unique RREF;
hull dimensions use stacked right-dual equations, independently of the formula.
F_9 is F_3[u]/(u^2+1). No source program is used.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import prod


class Field:
    def __init__(self, p, d):
        assert p >= 3 and all(p % t for t in range(2, p))
        assert d == 1 or (p, d) == (3, 2)
        self.p, self.d, self.q = p, d, p**d
        q = self.q
        if d == 1:
            self.add = [[(a+b) % p for b in range(q)] for a in range(q)]
            self.mul = [[a*b % p for b in range(q)] for a in range(q)]
        else:
            self.add = [[(a % p+b % p) % p+p*((a//p+b//p) % p)
                         for b in range(q)] for a in range(q)]
            self.mul = [[((a % p)*(b % p)-(a//p)*(b//p)) % p
                         +p*(((a % p)*(b//p)+(a//p)*(b % p)) % p)
                         for b in range(q)] for a in range(q)]
        self.neg = [next(b for b in range(q) if self.add[a][b] == 0)
                    for a in range(q)]
        self.inv = [0]+[next(b for b in range(1, q) if self.mul[a][b] == 1)
                        for a in range(1, q)]
        self.frob = [[self.power(a, p**i) for a in range(q)] for i in range(d)]
        for a in range(q):
            assert self.power(a, q) == a
            for b in range(q):
                for i in range(d):
                    f = self.frob[i]
                    assert f[self.add[a][b]] == self.add[f[a]][f[b]]
                    assert f[self.mul[a][b]] == self.mul[f[a]][f[b]]

    def power(self, a, exponent):
        result = 1
        while exponent:
            if exponent % 2:
                result = self.mul[result][a]
            a = self.mul[a][a]
            exponent //= 2
        return result

    def dot(self, x, y):
        result = 0
        for a, b in zip(x, y):
            result = self.add[result][self.mul[a][b]]
        return result

    def twist(self, x, i):
        return tuple(self.frob[i % self.d][a] for a in x)

    def rref(self, rows, n):
        rows = [list(x) for x in rows]
        pivot = 0
        for col in range(n):
            j = next((j for j in range(pivot, len(rows)) if rows[j][col]), None)
            if j is None:
                continue
            rows[pivot], rows[j] = rows[j], rows[pivot]
            a = self.inv[rows[pivot][col]]
            rows[pivot] = [self.mul[a][v] for v in rows[pivot]]
            for j in range(len(rows)):
                if j != pivot and rows[j][col]:
                    a = self.neg[rows[j][col]]
                    rows[j] = [self.add[x][self.mul[a][y]]
                               for x, y in zip(rows[j], rows[pivot])]
            pivot += 1
            if pivot == len(rows):
                break
        return tuple(tuple(x) for x in rows[:pivot])

    def span(self, rows, n):
        for coefficients in product(range(self.q), repeat=len(rows)):
            v = [0]*n
            for a, row in zip(coefficients, rows):
                for j in range(n):
                    v[j] = self.add[v[j]][self.mul[a][row[j]]]
            yield tuple(v)


def subspaces(F, n, k):
    """Every k-subspace, with no generator-matrix multiplicity."""
    for pivots in combinations(range(n), k):
        free = [(i, j) for i, p in enumerate(pivots)
                for j in range(p+1, n) if j not in pivots]
        for values in product(range(F.q), repeat=len(free)):
            rows = [[int(j == p) for j in range(n)] for p in pivots]
            for (i, j), x in zip(free, values):
                rows[i][j] = x
            yield tuple(tuple(row) for row in rows)


@lru_cache(None)
def gaussian(q, a, b):
    if not 0 <= b <= a:
        return 0
    result = Fraction(prod(q**(a-i)-1 for i in range(b)),
                      prod(q**i-1 for i in range(1, b+1)))
    assert result.denominator == 1
    return result.numerator


def mobius(p, u):
    return (-1)**u*p**(u*(u-1)//2)


@lru_cache(None)
def isotropic(p, n, r):
    if not 0 <= r <= n//2:
        return 0
    m = n//2
    if n % 2:
        numerator = prod(p**(2*(m-u))-1 for u in range(r))
    else:
        eps = 1 if pow((-1)**m % p, (p-1)//2, p) == 1 else -1
        numerator = prod((p**(m-u)-eps)*(p**(m-u-1)+eps)
                         for u in range(r))
    result = Fraction(numerator, prod(p**i-1 for i in range(1, r+1)))
    assert result.denominator == 1 and result >= 0
    return result.numerator


@lru_cache(None)
def exact_closure(p, d, r, j):
    return sum(mobius(p, r-t)*gaussian(p, r, t)*gaussian(p**d, t, j)
               for t in range(r+1))


def moment(p, d, n, k, j):
    return sum(isotropic(p, n, r)*exact_closure(p, d, r, j)
               *gaussian(p**d, n-r-j, k-j)
               for r in range(j, min(d*j, n//2)+1))


def mass(p, d, n, k, h):
    q = p**d
    return sum(mobius(q, j-h)*gaussian(q, j, h)*moment(p, d, n, k, j)
               for j in range(h, k+1))


def hull_dimension(F, basis):
    k = len(basis)
    # For y=sum_l a_l b_l, apply phi^(-i) to B(x,phi^i(y))=0.
    equations = [[F.frob[(-i) % F.d][F.dot(x, F.twist(y, i))]
                  for y in basis] for i in range(F.d) for x in basis]
    return k-len(F.rref(equations, k))


def literal_hull_size(F, basis, n):
    # Direct set definition, without linearizing the semilinear equations.
    return sum(all(F.dot(x, F.twist(y, i)) == 0
                   for i in range(F.d) for x in basis)
               for y in F.span(basis, n))


def fidelity(F):
    # Source Lemma 45, all ordered vector pairs in F_q^2. The final comparison
    # tests Euclidean orthogonality of the spans of all Frobenius conjugates.
    vectors = list(product(range(F.q), repeat=2))
    checked = 0
    for x in vectors:
        for y in vectors:
            total = all(F.dot(x, F.twist(y, i)) == 0 for i in range(F.d))
            closure_euclidean = all(F.dot(F.twist(x, a), F.twist(y, b)) == 0
                                    for a in range(F.d) for b in range(F.d))
            closure_total = all(F.dot(F.twist(x, a), F.twist(y, b+i)) == 0
                                for a in range(F.d) for b in range(F.d)
                                for i in range(F.d))
            assert total == closure_euclidean == closure_total, (x, y)
            checked += 1
    print(f'FIDELITY Lemma45 p={F.p} d={F.d} n=2 ordered_pairs={checked}')


def check_family(p, d, max_n):
    F = Field(p, d)
    fidelity(F)
    cases = codes = 0
    for n in range(max_n+1):
        for k in range(n+1):
            observed = [0]*(k+1)
            closure_observed = [0]*(min(d*k, n)+1)
            for index, basis in enumerate(subspaces(F, n, k)):
                h = hull_dimension(F, basis)
                observed[h] += 1
                closure = F.rref([F.twist(x, i) for x in basis
                                  for i in range(d)], n)
                closure_observed[len(closure)] += 1
                # Full source set-definition audit on n<=2; deterministic
                # samples for n=3,4, including the full ambient space.
                if n <= 2 or index < 12:
                    assert literal_hull_size(F, basis, n) == F.q**h
                # Source Lemma45 at subspace level, checked for every code.
                assert (h == k) == all(F.dot(x, y) == 0
                                      for x in closure for y in closure)
            expected = [mass(p, d, n, k, h) for h in range(k+1)]
            assert observed == expected, (p, d, n, k, observed, expected)
            assert sum(observed) == gaussian(F.q, n, k)
            # Separate closure-count and hull-incidence identities.
            for r, actual in enumerate(closure_observed):
                assert actual == gaussian(p, n, r)*exact_closure(p, d, r, k)
            for j in range(k+1):
                assert moment(p, d, n, k, j) == sum(
                    gaussian(F.q, ell, j)*observed[ell] for ell in range(j, k+1))
            if d == 1:
                assert observed[k] == isotropic(p, n, k)
            print(f'MASS p={p} d={d} n={n} k={k} codes={sum(observed)} M={observed}')
            cases += len(observed)
            codes += sum(observed)
    return cases, codes


def main():
    cases = codes = 0
    for p, d, max_n in [(3, 1, 4), (3, 2, 4), (5, 1, 3)]:
        c, s = check_family(p, d, max_n)
        cases += c
        codes += s
    print(f'SUMMARY mass_parameter_tuples={cases} enumerated_subspaces={codes}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
