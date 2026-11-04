#!/usr/bin/env python3
"""Exact exhaustive small-case factorization check for arXiv:2610.23684v2.

All divisors of (x^a+x^-a)^n have profiles (x^a+x^-a)^k. The README
explains the finite reduction, including arbitrary real exponents. We enumerate
all admissible fixed-point counts, test all proper binary factorizations to
find atoms, then enumerate every unordered multiset of those atoms.
"""
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import product
from math import comb
import sympy as sp


def multiply(p, q):
    out = Counter()
    for u, c in p.items():
        for v, d in q.items():
            out[u + v] += c * d
    return dict(out)


def power(p, n):
    out = {Fraction(0): 1}
    for _ in range(n):
        out = multiply(out, p)
    return out


def flags(k):
    constant = comb(k, k // 2) if k % 2 == 0 else 0
    return tuple(f for f in range(constant + 1) if f % 2 == constant % 2)


def all_atoms(n):
    # Every split 1 <= i < k has strictly smaller nonunit profile degree.
    atoms = []
    for k in range(1, n + 1):
        for f in flags(k):
            reducible = any(g * h == f
                            for i in range(1, k)
                            for g in flags(i)
                            for h in flags(k - i))
            if not reducible:
                atoms.append((k, f))
    return tuple(atoms)


def census(n, target_f=0):
    atoms = all_atoms(n)

    @cache
    def visit(remaining, start, f_product):
        if remaining == 0:
            return ((0, 1),) if f_product == target_f else ()
        counts = Counter()
        for idx in range(start, len(atoms)):
            k, f = atoms[idx]
            if k > remaining:
                break
            for ell, count in visit(remaining - k, idx, f_product * f):
                counts[ell + 1] += count
        return tuple(sorted(counts.items()))

    return dict(visit(n, 0, 1)), atoms


def verify_profiles(a, n):
    base = {a: 1, -a: 1}
    profiles = [power(base, k) for k in range(n + 1)]
    for k, p in enumerate(profiles):
        # Independent expansion by all choices of the two orbit sites.
        brute = Counter(sum(xs, Fraction(0)) for xs in product((a, -a), repeat=k))
        assert p == dict(brute)
        assert p == {(2 * i - k) * a: comb(k, i) for i in range(k + 1)}
        assert all(p.get(-u, 0) == c for u, c in p.items())
        assert sum(p.values()) == 2 ** k
        assert p.get(Fraction(0), 0) == (comb(k, k // 2) if k % 2 == 0 else 0)
        assert flags(k)
    for k in range(n + 1):
        for i in range(k + 1):
            assert multiply(profiles[i], profiles[k - i]) == profiles[k]
    return profiles


def main():
    z = sp.Symbol('z')
    for n in range(1, 12):
        coefficient, irreducible_factors = sp.factor_list(sp.Poly((1 + z) ** n, z))
        assert coefficient == 1
        assert irreducible_factors == [(sp.Poly(z + 1, z), n)]
    print('NORMALIZED_DIVISORS: exact factor_list((1+z)^n), n=1..11: only (z+1)^k')
    print('METHOD: all admissible fixed-point flags; all binary splits; all unordered atom multisets')
    for t in (Fraction(2), Fraction(3)):
        a = t / 2
        profiles = verify_profiles(a, 11)
        for j in range(1, 6):
            n = 2 * j + 1
            counts, atoms = census(n)
            lengths = sorted(counts)
            assert lengths == list(range(2, n + 1))
            # Fidelity: source Lemma 4.17 / Theorem 4.18 and Remark 4.20.
            assert (1, 0) in atoms
            for h in range(1, j + 1):
                assert 2 in flags(2 * h)
                assert (2 * h, 2) in atoms
            assert multiply(profiles[2 * j], profiles[1]) == profiles[n]
            assert counts[2] > 0 and counts[n] > 0
            if j <= 3:
                print(f'FIDELITY Remark 4.20: t={t}, j={j}, lengths={lengths}')
            # Check every proposed witness, not just the target length set.
            witnesses = [[(1, 0)] * n]
            witnesses += [[(2 * (j + 1 - r), 2)] + [(1, 0)] * (2 * r - 1)
                          for r in range(1, j + 1)]
            witnesses += [[(2 * (j - r), 2), (2, 2)] + [(1, 0)] * (2 * r - 1)
                          for r in range(1, j)]
            assert sorted(map(len, witnesses)) == list(range(2, n + 1))
            for factors in witnesses:
                p, f = {Fraction(0): 1}, 1
                for k, flag in factors:
                    assert (k, flag) in atoms
                    p, f = multiply(p, profiles[k]), f * flag
                assert p == profiles[n] and f == 0
            print(f'EXHAUSTIVE t={t}, a={a}, j={j}, sites={2 ** n}, '
                  f'atom_types={len(atoms)}, factorizations_by_length={counts}, lengths={lengths}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
