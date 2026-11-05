#!/usr/bin/env python3
"""Exact finite checks for Ashtab--Fernandez-Breton, arXiv:2607.26522, Q5.2(1).

Periodic-algebra models test the proposed support/stabilizer construction.
Finite windows test its averaging estimates without approximating an ultralimit.
These checks do not prove existence on the full powerset of the integers.
"""
from fractions import Fraction as Q
from itertools import product


def subsets(n):
    return range(1 << n)


def members(mask, n):
    return {i for i in range(n) if mask & (1 << i)}


def source_parity_check():
    # Section 4, immediately after Definition 4.1: every invariant probability
    # assigns equal mass to the even and odd integers, hence mass 1/2 to each.
    for r in range(1, 33):
        m = 2 * r
        even = {i for i in range(m) if i % 2 == 0}
        odd = {i for i in range(m) if i % 2 == 1}
        base = lambda a: Q(len(a), m)
        assert even.isdisjoint(odd) and even | odd == set(range(m))
        assert {(x + 1) % m for x in even} == odd
        assert base(even) == base(odd) == Q(1, 2)
        assert base(even) + base(odd) == 1
    print('FIDELITY: source parity identity m(2Z) = m(2Z+1) = 1/2; 32 periodic algebras')


def periodic_family_check():
    models = sets = shift_checks = 0
    for q in range(2, 7):
        for r in range(2, 7):
            m = q * r
            if m > 12:
                continue
            models += 1
            universe = set(range(m))
            h = {x for x in universe if x % q == 0}
            mu = lambda a: Q(len(a & h), len(h))
            translate = lambda a, x: {(z - x) % m for z in a}
            assert mu(h) == mu(universe) == 1 and mu(set()) == 0
            assert mu(translate(h, 1)) == 0
            # The r residue classes in H all have mass 1/r, so non-binary.
            assert mu({0}) == Q(1, r) and mu({0}) not in (0, 1)
            if r % 2 == 0:
                double_h = {x for x in universe if x % (2 * q) == 0}
                assert mu(double_h) == Q(1, 2)
            for mask in subsets(m):
                a = members(mask, m)
                sets += 1
                value = mu(a)
                assert 0 <= value <= 1
                # Independently compute the same measure by conditioning the
                # uniform probability on the finite quotient on H.
                assert value == q * Q(len(a & h), m)
                good = set()
                for x in universe:
                    shift_checks += 1
                    shifted = mu(translate(a, x))
                    if x in h:
                        assert shifted == value
                    if shifted == value:
                        good.add(x)
                assert h <= good and mu(good) == 1
    print(f'FAMILY: q=2..6, r=2..6, qr<=12; {models} models, {sets} subsets, {shift_checks} shifts')
    print('FAMILY: support mass 1, subgroup invariance, good-shift mass 1, both separations all hold')


def finite_additivity_check():
    checks = 0
    for n in range(1, 8):
        # Every disjoint pair in an n-point sample is represented once.
        for tags in product(range(3), repeat=n):
            a = {i for i, tag in enumerate(tags) if tag == 1}
            b = {i for i, tag in enumerate(tags) if tag == 2}
            assert a.isdisjoint(b)
            assert Q(len(a | b), n) == Q(len(a), n) + Q(len(b), n)
            checks += 1
    print(f'AVERAGES: exact finite additivity on {checks} disjoint pairs')


def window_check():
    checks = 0
    for n in range(1, 7):
        initial = set(range(1, n + 1))
        for k in range(-3, 4):
            shifted = set(range(1 + k, n + 1 + k))
            joint = sorted(initial | shifted)
            assert len(initial ^ shifted) == 2 * min(n, abs(k))
            for mask in subsets(len(joint)):
                a = {z for j, z in enumerate(joint) if mask & (1 << j)}
                gap = abs(Q(len(a & initial), n) - Q(len(a & shifted), n))
                assert gap <= Q(len(initial ^ shifted), n) <= Q(2 * abs(k), n)
                # A stronger bound, from matching the equal-sized windows.
                assert gap <= Q(min(n, abs(k)), n)
                checks += 1
    for q in range(2, 7):
        for n in range(1, 65):
            sample = {q * j for j in range(1, n + 1)}
            assert len(sample) == n
            assert all(z % q == 0 for z in sample)
            assert all((z + 1) % q != 0 for z in sample)
            double = Q(sum(z % (2 * q) == 0 for z in sample), n)
            assert double == Q(n // 2, n)
            assert abs(double - Q(1, 2)) <= Q(1, n)
    print(f'WINDOWS: {checks} exhaustive subset/shift checks; N=1..6, k=-3..3')
    print('WINDOWS: q=2..6, N=1..64; H mass 1, H-1 mass 0, 2H average floor(N/2)/N')
    print('SCOPE: finite exact checks; full-powerset ultrafilter-limit existence is not a computational result')


def main():
    source_parity_check()
    periodic_family_check()
    finite_additivity_check()
    window_check()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
