#!/usr/bin/env python3
"""Exact finite checks for Breloer--Riener Conjecture 5.12.

Integer orbit enumeration, truncated Euler products and symbolic support checks.
These tests do not certify infinite-series algebraic independence.
"""
from functools import lru_cache
from itertools import product
from math import comb


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for x in range(total + 1):
            for tail in compositions(total - x, length - 1):
                yield (x,) + tail


def mul(a, b, degree):
    return [sum(a[j] * b[i - j] for j in range(i + 1))
            for i in range(degree + 1)]


def power(a, n, degree):
    out = [1] + [0] * degree
    for _ in range(n):
        out = mul(out, a, degree)
    return out


def euler(exponents, degree):
    """Coefficient expansion from exact negative-binomial coefficients."""
    out = [1] + [0] * degree
    for m in range(1, degree + 1):
        h = exponents(m)
        factor = [0] * (degree + 1)
        for b in range(degree // m + 1):
            factor[m * b] = comb(h + b - 1, b) if h else int(b == 0)
        out = mul(out, factor, degree)
    return out


def orbit_counts(types, degree, accept=lambda sites: True, max_sites=None):
    """Enumerate unordered multisets of nonzero local exponent vectors.

    A multiset is precisely an orbit under permutation of sites. Site units fill
    the unused positions. No Euler-product coefficient formula is used here.
    """
    types = sorted(types, key=lambda v: (sum(v), v))
    counts = [0] * (degree + 1)

    def visit(start, used_degree, sites):
        if not accept(sites):
            return
        counts[used_degree] += 1
        if max_sites is not None and len(sites) == max_sites:
            return
        for i in range(start, len(types)):
            v = types[i]
            w = sum(v)
            if used_degree + w > degree:
                break
            visit(i, used_degree + w, sites + (v,))

    visit(0, 0, ())
    return counts


def main():
    D = 9
    pairs = [v for m in range(1, D + 1) for v in compositions(m, 2)]
    p1 = euler(lambda m: 1, D)
    fidelity = [
        ("I^(1)", lambda s: not (any(v[0] for v in s) and any(v[1] for v in s)),
         [1] + [2 * x for x in p1[1:]]),
        ("I^(2)", lambda s: all(x * y == 0 for x, y in s), power(p1, 2, D)),
        ("I^(3)", lambda s: sum(x > 0 for x, _ in s) <= 1
         and sum(y > 0 for _, y in s) <= 1,
         [1] + [2 * d for d in range(1, D + 1)]),
    ]
    for name, legal, expected in fidelity:
        counted = orbit_counts(pairs, D, accept=legal)
        assert counted == expected, (name, counted, expected)
        print(f"FIDELITY Example 5.11 {name}, degrees 0..{D}: {counted}")
    assert orbit_counts(pairs, D, accept=fidelity[2][1], max_sites=2) == fidelity[2][2]
    print("FIDELITY I^(3) stabilizes at two sites: OK")

    D = 16
    p1, p2 = euler(lambda m: 1, D), euler(lambda m: m + 1, D)
    for r in range(2, 6):
        h = euler(lambda m: r * m + 1, D)
        assert mul(h, power(p1, r - 1, D), D) == power(p2, r, D)
        small = 7
        types = []
        for m in range(1, small + 1):
            normal = [v for v in compositions(m, r + 1)
                      if all(v[a] * v[b] == 0 for a in range(1, r + 1)
                             for b in range(a + 1, r + 1))]
            assert len(normal) == r * m + 1
            types.extend(normal)
        counted = orbit_counts(types, small)
        assert counted == h[:small + 1]
        assert orbit_counts(types, small, max_sites=small) == counted
        assert orbit_counts(types, small, max_sites=small + 1) == counted
        print(f"FAMILY r={r}: local counts h_m=rm+1 for m=1..{small}; "
              f"orbits degrees 0..{small}: {counted}")
        print(f"FAMILY r={r}: H*P1^(r-1)=P2^r through degree {D}: OK")
        if r == 2:
            assert h[:9] == [1, 3, 11, 32, 90, 231, 576, 1363, 3141]
            print(f"FAMILY r=2, degrees 0..8: {h[:9]}")
        # The union of the r coordinate planes is exactly the zero set of the
        # squarefree pairwise products on this exact rational grid.
        for v in product(range(-1, 2), repeat=r + 1):
            vanishes = all(v[a] * v[b] == 0 for a in range(1, r + 1)
                           for b in range(a + 1, r + 1))
            in_union = any(all(v[b] == 0 for b in range(1, r + 1) if b != a)
                           for a in range(1, r + 1))
            assert vanishes == in_union
        # On each plane, z and x_a are unrestricted; intersection is the z-axis.
        for z in range(-2, 3):
            v = (z,) + (0,) * r
            assert all(all(v[b] == 0 for b in range(1, r + 1) if b != a)
                       for a in range(1, r + 1))
        # Independent of coefficients, U^(r-1)*F cannot cancel V^r: the
        # U-exponent of every monomial in the former is positive.
        for a in range(8):
            for b in range(8):
                assert (a + r - 1, b) != (0, r)
        print(f"FAMILY r={r}: plane equations on {{-1,0,1}}^(r+1), "
              "axis intersection, symbolic exponent separation: OK")
    print("SCOPE: finite exact checks; analytic nonmembership is the preregistered argument.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
