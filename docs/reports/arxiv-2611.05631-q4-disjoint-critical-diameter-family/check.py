#!/usr/bin/env python3
"""Exact finite checks for the Lopez-Campos Question 4 candidate.

The source fidelity check is Proposition 1, every k-chromatic graph contains
an induced vertex-k-critical subgraph, on every labelled graph on <=5 vertices.
The family checks are finite diagnostics, not a proof of the infinite family.
Only integer and Fraction arithmetic is used. No external packages are needed.
"""
from fractions import Fraction
from functools import cache
from itertools import combinations, combinations_with_replacement
from math import comb, prod, isqrt


def graph_checks():
    graph_count = packing_count = 0
    for n in range(6):
        edges = list(combinations(range(n), 2))
        for code in range(1 << len(edges)):
            adj = [0] * n
            for bit, (u, v) in enumerate(edges):
                if code >> bit & 1:
                    adj[u] |= 1 << v
                    adj[v] |= 1 << u

            def independent(mask):
                return all(not (mask & adj[u]) for u in range(n) if mask >> u & 1)

            @cache
            def chi(mask):
                if not mask:
                    return 0
                least = mask & -mask
                best, sub = n + 1, mask
                while sub:
                    if sub & least and independent(sub):
                        best = min(best, 1 + chi(mask ^ sub))
                    sub = (sub - 1) & mask
                return best

            def extract(mask, k):
                assert chi(mask) >= k > 0
                for u in range(n):
                    smaller = mask & ~(1 << u)
                    if chi(smaller) >= k:
                        mask = smaller
                assert chi(mask) == k
                assert all(chi(mask & ~(1 << u)) < k
                           for u in range(n) if mask >> u & 1)
                return mask

            whole = (1 << n) - 1
            c = chi(whole)
            if c:
                extract(whole, c)  # Source Proposition 1, with induced criticality.
            for k in range(1, c + 1):
                for t in range(1, c // k + 1):
                    remainder, used = whole, 0
                    for i in range(t):
                        assert chi(remainder) >= (t - i) * k
                        block = extract(remainder, k)
                        assert not (block & used)
                        assert chi(remainder) <= k + chi(remainder ^ block)
                        used |= block
                        remainder ^= block
                    packing_count += 1
            graph_count += 1
    print(f"FIDELITY source Proposition 1: all {graph_count} labelled graphs on 0..5 vertices")
    print(f"PACKING chi(G)>=t*k: {packing_count} admissible (G,k,t) cases, all passed")


def multilinear_polynomial(a, p):
    # Work in F_p[X]/(X_i^2-X_i); monomials are bit masks.
    terms = {0: 1}
    for r in range(1, p):
        updated = {}
        for monomial, coefficient in terms.items():
            updated[monomial] = (updated.get(monomial, 0) - r * coefficient) % p
            for i in a:
                target = monomial | (1 << i)
                updated[target] = (updated.get(target, 0) + coefficient) % p
        terms = {m: c for m, c in updated.items() if c}
    assert all(m.bit_count() <= p - 1 for m in terms)
    return terms


def geometry_checks(p):
    m = 4 * p
    sets = [frozenset((0,) + rest) for rest in combinations(range(1, m), 2 * p - 1)]
    masks = [sum(1 << i for i in a) for a in sets]
    pairs = list(combinations(range(m), 2))
    numerators = []
    for a in sets:
        x = [1 if i in a else -1 for i in range(m)]
        numerators.append(tuple(x[i] * x[j] for i, j in pairs))
    assert len(set(numerators)) == len(sets) == comb(m, 2 * p) // 2
    polynomials = [multilinear_polynomial(a, p) for a in sets] if p > 1 else []
    max_sq, checked, edge_count = Fraction(0), 0, 0
    for ai, bi in combinations_with_replacement(range(len(sets)), 2):
        a, b = sets[ai], sets[bi]
        k = len(a & b)
        lhs = Fraction(sum((u - v) ** 2 for u, v in zip(numerators[ai], numerators[bi])), m * m)
        rhs = 1 - Fraction((k - p) ** 2, p * p)
        assert lhs == rhs
        assert 0 <= lhs <= 1
        assert (ai != bi and lhs == 1) == (k == p)
        max_sq = max(max_sq, lhs)
        edge_count += int(ai != bi and lhs == 1)
        if p > 1:
            value = prod(k - r for r in range(1, p)) % p
            ml_value = sum(c for monomial, c in polynomials[ai].items()
                           if monomial & masks[bi] == monomial) % p
            assert value == ml_value
            if ai == bi:
                assert value != 0
            elif k != p:
                assert value == 0
        checked += 1
    assert max_sq == 1
    print(f"GEOMETRY p={p}: vertices={len(sets)}, dimension={comb(m,2)}, pairs_with_diagonal={checked}, diameter_squared=1, edges={edge_count}")
    if p > 1:
        print(f"POLYNOMIAL p={p}: all tested evaluations and degree<={p-1} multilinearizations passed")


def multiply(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def h_shift(base):
    # h(p)=(4p+1)(4p)(4p-1); p=base+q.
    return multiply(multiply([4 * base + 1, 4], [4 * base, 4]), [4 * base - 1, 4])


def bound_checks():
    primes = [p for p in range(37, 101)
              if all(p % q for q in range(2, isqrt(p) + 1))]
    for p in primes:
        d = comb(4 * p, 2)
        vertices = comb(4 * p, 2 * p) // 2
        dimension_bound = sum(comb(4 * p, j) for j in range(p))
        assert comb(4 * p, 2 * p) * (4 * p + 1) >= 16 ** p
        assert 3 * dimension_bound * 27 ** p <= 256 ** p
        lower = Fraction(3 * 27 ** p, 2 * (4 * p + 1) * 16 ** p)
        assert Fraction(vertices, dimension_bound) >= lower > 2 * d
        supplier_n = 4 * p - 2
        supplier_d = comb(supplier_n, 2)
        supplier_vertices = 2 ** (supplier_n - 2)
        supplier_bound = sum(comb(supplier_n - 1, j) for j in range(p - 1))
        assert supplier_d <= d
        assert 64 * supplier_bound * 27 ** p <= 3 * 256 ** p
        assert Fraction(supplier_vertices, supplier_bound) >= Fraction(4, 3) * Fraction(27, 16) ** p
        assert Fraction(supplier_vertices, supplier_bound) > 2 * d
        if p == 37:
            print(f"REUSE FormalBook borsuk_theorem p=37: source_dimension={supplier_d}, target_dimension={d}, ceil(partition_bound)={(supplier_vertices+supplier_bound-1)//supplier_bound}")
            print(f"BOUND p=37: d={d}, 2*d={2*d}")
            print(f"BOUND p=37: L=3*(27/16)^37/(2*149)={lower}")
            print(f"BOUND p=37: L>2*d, cross-multiplied margin={3*27**p-4*(4*p+1)*d*16**p}")
            print(f"BOUND p=37: ceil(vertices/dimension_bound)={(vertices+dimension_bound-1)//dimension_bound}")
    # This positive-coefficient identity is the induction step for all p>=37.
    coefficients = [27 * a - 16 * b for a, b in zip(h_shift(37), h_shift(38))]
    assert all(c > 0 for c in coefficients)
    for q in range(20):
        value = sum(c * q ** i for i, c in enumerate(coefficients))
        p = q + 37
        assert value == 27 * (4*p+1)*(4*p)*(4*p-1) - 16 * (4*p+5)*(4*p+4)*(4*p+3)
    print(f"BOUND all p>=37 induction polynomial, coefficients in ascending powers of (p-37): {coefficients}")
    print(f"BOUND exact binomial estimates additionally checked for primes {primes}")


if __name__ == '__main__':
    graph_checks()
    for small_p in (1, 2, 3):
        geometry_checks(small_p)
    bound_checks()
    print('ALL_CHECKS_OK')
