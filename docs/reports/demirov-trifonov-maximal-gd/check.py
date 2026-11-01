#!/usr/bin/env python3
"""Exact finite checks for Demirov--Trifonov Problem 7, Examples 2 and 3.

Uses Python's integers only. These checks do not prove global maximality
against arbitrary infinite, potentially nonperiodic integer extensions.
"""
from itertools import combinations
from math import comb


def primes_through(bound):
    return [p for p in range(2, bound + 1)
            if all(p % q for q in range(2, p) if q * q <= p)]


def trim(poly):
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def remainder(a, b, p):
    a = trim([x % p for x in a])
    b = trim([x % p for x in b])
    assert b
    inv_lead = pow(b[-1], -1, p)
    while len(a) >= len(b):
        shift = len(a) - len(b)
        factor = a[-1] * inv_lead % p
        for j, coeff in enumerate(b):
            a[shift + j] = (a[shift + j] - factor * coeff) % p
        trim(a)
    return a


def gcd_degree(d, p):
    # X^d - 1 and (1-X)^d - 1 over F_p, in ascending coefficient order.
    f = [(-1) % p] + [0] * (d - 1) + [1]
    g = [((-1) ** i * comb(d, i)) % p for i in range(d + 1)]
    g[0] = (g[0] - 1) % p
    g = trim(g)
    while g:
        f, g = g, remainder(f, g, p)
    assert f
    return len(f) - 1


def sum_free(residues, modulus):
    return all((a + b) % modulus not in residues
               for a in residues for b in residues)


def multiplicatively_closed(residues, modulus):
    return all(a * b % modulus in residues
               for a in residues for b in residues)


def divisor_roots(p):
    return {d: frozenset(a for a in range(1, p) if pow(a, d, p) == 1)
            for d in range(1, p) if (p - 1) % d == 0}


def source_fidelity():
    # Theorem 5(i): V_k has a proper GD extension by the +/-1 classes, k>=4.
    # The entire finite quotient is checked, not just an initial integer segment.
    for k in range(4, 101):
        v = {1}
        extension = {1, k - 1}
        assert v < extension
        assert multiplicatively_closed(extension, k)
        assert sum_free(extension, k)
    # Source Examples 2 and 3, with source indexing n>=1 converted to j>=0.
    for k in range(2, 51):
        for i in range(12):
            for j in range(12):
                a, b = 1 + k * i, 1 + k * j
                assert (a * b - 1) % k == 0
                assert (a + b - 1) % k != 0
    for k in range(3, 31):
        powers = {k ** j for j in range(25)}
        for i in range(12):
            for j in range(12):
                assert k ** i * k ** j == k ** (i + j)
                assert k ** i + k ** j not in powers
    print('FIDELITY: Theorem 5(i) strict +/-1 residue extensions, 4<=k<=100: OK')
    print('FIDELITY: Example 2, 2<=k<=50; Example 3, 3<=k<=30, exponents 0..11: OK')


def check_construction():
    selected = {}
    tested_orders = 0
    for p in primes_through(199):
        roots = divisor_roots(p)
        good = []
        for d, h in roots.items():
            assert len(h) == d
            assert multiplicatively_closed(h, p)
            direct = sum_free(h, p)
            normalized = all((1 - a) % p not in h for a in h)
            polynomial = gcd_degree(d, p) == 0
            assert direct == normalized == polynomial, (p, d)
            tested_orders += 1
            if direct:
                good.append(d)
        assert 1 in good
        d = max(good)
        h = roots[d]
        assert 1 in h and 0 not in h
        assert all(not (h < roots[e]) for e in good)
        ordered = sorted(h)
        formula = [p * (j // d) + ordered[j % d] for j in range(6 * d)]
        scanned = [a for a in range(1, 6 * p) if pow(a, d, p) == 1]
        assert formula == scanned
        assert all(a < b for a, b in zip(formula, formula[1:]))
        # Exhaustive one-element generated subgroup extension obstruction.
        for residue in set(range(1, p)) - h:
            extension = set(h) | {residue}
            while True:
                enlarged = extension | {a * b % p for a in extension for b in extension}
                if enlarged == extension:
                    break
                extension = enlarged
            assert not sum_free(extension, p)
        selected[p] = (d, h)
        print(f'CONSTRUCTION p={p}: d={d}; roots={ordered}; gcd/residue/maximal-subgroup/term-formula OK')
    print(f'EXACT ALGEBRA: {tested_orders} divisor orders for all {len(selected)} primes <=199: OK')
    return selected


def check_finite_classification():
    checked = 0
    for p in primes_through(13):
        known = set(divisor_roots(p).values())
        closed = set()
        for mask in range(1 << (p - 2)):
            h = frozenset({1} | {r for r in range(2, p) if mask & (1 << (r - 2))})
            checked += 1
            if multiplicatively_closed(h, p):
                assert h in known
                closed.add(h)
        assert closed == known
        maximal_good = [h for h in closed if sum_free(h, p)
                        and not any(h < k and sum_free(k, p) for k in closed)]
        assert maximal_good
    print(f'FINITE CLASSIFICATION: all {checked} residue subsets containing 1, primes <=13: OK')


def check_lifting_identity():
    checked = 0
    for p in primes_through(43):
        for a in range(1, p):
            for b in range(1, p):
                c0 = (a + b) % p
                if not c0:
                    continue
                # Different positive representatives, to test the integer lift itself.
                for a1, b1, c in [(a, b, c0), (a + p, b + 2 * p, c0 + 3 * p)]:
                    u = pow(b1 * c, -1, p)
                    v = pow(a1 * c, -1, p)
                    assert u > 0 and v > 0
                    assert b1 * c * u % p == 1
                    assert a1 * c * v % p == 1
                    assert a1 * b1 * (u + v) % p == 1
                    left = a1 * (b1 * c * u)
                    right = b1 * (a1 * c * v)
                    total = c * (a1 * b1 * (u + v))
                    assert left + right == total
                    checked += 1
    print(f'INTEGER LIFT: {checked} exact congruence/addition certificates, primes <=43: OK')


def check_family_inclusions(selected):
    for k in range(2, 151):
        p = min(q for q in selected if k % q == 0)
        d, _ = selected[p]
        for j in range(101):
            assert pow(1 + k * j, d, p) == 1
    for k in range(3, 151):
        p = min(q for q in selected if (k - 1) % q == 0)
        d, _ = selected[p]
        for j in range(101):
            assert pow(pow(k, j), d, p) == 1
    print('FAMILY INCLUSIONS: V_k, 2<=k<=150; W_k, 3<=k<=150; j=0..100: OK')


def main():
    source_fidelity()
    selected = check_construction()
    check_finite_classification()
    check_lifting_identity()
    check_family_inclusions(selected)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
