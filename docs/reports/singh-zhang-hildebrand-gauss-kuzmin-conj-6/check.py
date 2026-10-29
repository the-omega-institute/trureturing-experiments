#!/usr/bin/env python3
"""Exact finite checks for Singh--Zhang--Hildebrand Conjecture 6.

Source: arXiv:2511.01992v1, Integers 25 (2025), A89.
Uses only integer arithmetic and Fraction. No logarithms are approximated:
the endpoint cross-ratio is the rational number whose base-two log is the measure.
Finite checks support the proposed argument; they do not prove the limiting claim.
"""
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, factorial


def convergent_matrix(word):
    """Direct product of source matrices [[0,1],[1,a]]."""
    p0, p, q0, q = 1, 0, 0, 1
    for a in word:
        p0, p, q0, q = p, p0 + a * p, q, q0 + a * q
    return p0, p, q0, q


def chi_matrix(word):
    _, p, q0, q = convergent_matrix(word)
    return (p + q) * (q0 + q)


def continuant(word):
    """Independent scalar recurrence, K(empty)=1, K(single a)=a."""
    previous, current = 0, 1
    for a in word:
        previous, current = current, a * current + previous
    return current


def chi_continuant(word):
    k = continuant(word)
    return (k + continuant(word[:-1])) * (k + continuant(word[1:]))


def add(p, q):
    out = p.copy()
    for e, c in q.items():
        out[e] = out.get(e, 0) + c
        if not out[e]:
            del out[e]
    return out


def mul(p, q):
    out = {}
    for e, c in p.items():
        for f, d in q.items():
            exponent = tuple(x + y for x, y in zip(e, f))
            out[exponent] = out.get(exponent, 0) + c * d
    return {e: c for e, c in out.items() if c}


def kpoly(indices, n):
    previous, current = {}, {(0,) * n: 1}
    for i in indices:
        e = tuple(int(j == i) for j in range(n))
        previous, current = current, add(mul({e: 1}, current), previous)
    return current


def chipoly(n):
    k = kpoly(range(n), n)
    left = kpoly(range(n - 1), n)
    right = kpoly(range(1, n), n)
    return mul(add(k, left), add(k, right))


def rename(poly, perm):
    out = {}
    for exponent, c in poly.items():
        e = [0] * len(perm)
        for i, j in enumerate(perm):
            e[j] = exponent[i]
        out[tuple(e)] = c
    return out


def evaluate(poly, values):
    total = 0
    for exponent, c in poly.items():
        term = c
        for a, power in zip(values, exponent):
            term *= a ** power
        total += term
    return total


def fidelity_and_numeric_bridge():
    assert chi_matrix((3, 1, 4)) == 552
    p0, p, q0, q = convergent_matrix((3, 1, 4))
    assert (Fraction(p, q), Fraction(p0 + p, q0 + q)) == (Fraction(5, 19), Fraction(6, 23))
    count = 0
    for n in range(1, 7):
        for word in product(range(1, 5), repeat=n):
            p0, p, q0, q = convergent_matrix(word)
            chi = chi_matrix(word)
            assert p0 * q - p * q0 == (-1) ** n
            assert chi == chi_continuant(word)
            assert chi == chi_matrix(word[::-1])
            u, v = Fraction(p, q), Fraction(p0 + p, q0 + q)
            ratio = max((1 + u) / (1 + v), (1 + v) / (1 + u))
            expected = Fraction(chi + 1, chi) if n % 2 == 0 else Fraction(chi, chi - 1)
            assert ratio == expected
            count += 1
    print(f"NUMERIC_BRIDGE: {count} words; lengths 1..6, digits 1..4; matrix, scalar continuant, reversal, determinant, endpoint ratio agree")
    print("FIDELITY: source illustration endpoints 5/19 and 6/23 agree; chi=552; exact ratio=552/551")
    assert Fraction(552, 551) != Fraction(551, 550)
    print("SOURCE_CORRECTION: the illustration rewrites the endpoint ratio as 551/550; that rewrite is false and is not used")


def source_tables_and_normalization():
    normalization_count = None
    for N, n, expected in [(10, 4, 10), (20, 4, 30), (10, 5, 8), (10, 6, 23)]:
        f = 0
        ordered_bad = 0
        for digits in combinations(range(1, N + 1), n):
            all_words = list(permutations(digits))
            values = {chi_matrix(w) for w in all_words}
            reps = [w for w in all_words if w < w[::-1]]
            rep_values = [chi_matrix(w) for w in reps]
            assert len(reps) == factorial(n) // 2
            assert set(rep_values) == values
            bad = len(values) < Fraction(factorial(n), 2)
            assert bad == (len(set(rep_values)) < len(rep_values))
            if bad:
                f += 1
                # Independently count exceptionality at each base ordering for N=10,n=4.
                if (N, n) == (10, 4):
                    for base in all_words:
                        orbit = {chi_matrix(tuple(base[i] for i in order))
                                 for order in permutations(range(n))}
                        ordered_bad += (len(orbit) < factorial(n) // 2)
        assert f == expected
        if (N, n) == (10, 4):
            assert ordered_bad == factorial(n) * f
            normalization_count = ordered_bad
        density = Fraction(f, comb(N, n))
        constant = comb(factorial(n) // 2, 2) * (2 * n - 1)
        falling = factorial(N) // factorial(N - n)
        assert density <= Fraction(constant * N ** (n - 1), falling)
        print(f"FIDELITY_TABLE: N={N}, n={n}, f={f}, delta={density}; source Tables 4/5/6 agree")
    print(f"NORMALIZATION: N=10,n=4, ordered_bad={normalization_count}, n!*f=240")


def polynomial_checks():
    for n in range(3, 8):
        poly = chipoly(n)
        k = kpoly(range(n), n)
        leading = (1,) * n
        assert k.get(leading) == 1
        assert all(sum(e) % 2 == n % 2 for e in k)
        expected_low = {}
        expected_high = {}
        for i in range(n - 1):
            e = tuple(0 if j in (i, i + 1) else 1 for j in range(n))
            expected_low[e] = 1
            e2 = tuple(1 if j in (i, i + 1) else 2 for j in range(n))
            expected_high[e2] = 2
        endpoints = tuple(1 if j in (0, n - 1) else 2 for j in range(n))
        expected_high[endpoints] = 1
        assert {e: c for e, c in k.items() if sum(e) == n - 2} == expected_low
        assert {e: c for e, c in poly.items() if sum(e) == 2 * n - 2} == expected_high
        assert {e: c for e, c in poly.items() if sum(e) == 2 * n} == {(2,) * n: 1}
        orbit = set()
        stabilizers = []
        for perm in permutations(range(n)):
            changed = rename(poly, perm)
            orbit.add(tuple(sorted(changed.items())))
            if changed == poly:
                stabilizers.append(perm)
            difference = add(changed, {e: -c for e, c in poly.items()})
            assert all(sum(e) <= 2 * n - 1 for e in difference)
        assert stabilizers == [tuple(range(n)), tuple(reversed(range(n)))]
        assert len(orbit) == factorial(n) // 2
        for word in product(range(1, 4), repeat=n):
            assert evaluate(poly, word) == chi_matrix(word)
        print(f"POLYNOMIAL: n={n}, monomials={len(poly)}, orbit={len(orbit)}, stabilizer=identity/reversal, coefficient pattern and comparison degree OK")


def main():
    fidelity_and_numeric_bridge()
    source_tables_and_normalization()
    polynomial_checks()
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
