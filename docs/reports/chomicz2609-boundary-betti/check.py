#!/usr/bin/env python3
"""Exact finite checks for Chomicz's secondary boundary-Betti question.

Source: arXiv:2609.04000v1, Example 2.9 and the paragraph after Lemma 3.13.
Only integer arithmetic is used. This finite check is not a proof of the
unbounded criterion. The source's Lemma 2.2 supplies completeness of the
three-pairing search for secondary semigroups.
"""
from functools import lru_cache
from itertools import combinations
from math import gcd, lcm

BOUND = 40


@lru_cache(maxsize=65536)
def factorizations(ds, n):
    """All nonnegative integer solutions, in the displayed generator order."""
    assert n >= 0 and all(d > 0 for d in ds)
    if len(ds) == 1:
        return ((n // ds[0],),) if n % ds[0] == 0 else ()
    return tuple((k,) + z for k in range(n // ds[0] + 1)
                 for z in factorizations(ds[1:], n - k * ds[0]))


def independent_factorizations(ds, n):
    """A separate four-generator enumeration with reversed loop order."""
    out = set()
    for z3 in range(n // ds[3] + 1):
        for z2 in range((n - z3 * ds[3]) // ds[2] + 1):
            rest = n - z3 * ds[3] - z2 * ds[2]
            for z1 in range(rest // ds[1] + 1):
                z0, rem = divmod(rest - z1 * ds[1], ds[0])
                if rem == 0:
                    out.add((z0, z1, z2, z3))
    return out


def graph_components(zs):
    """BFS in the graph joining distinct factorizations with shared support."""
    unseen = set(zs)
    out = set()
    while unseen:
        seed = unseen.pop()
        component = {seed}
        todo = [seed]
        while todo:
            z = todo.pop()
            neighbors = {w for w in unseen if any(a > 0 and b > 0
                                                  for a, b in zip(z, w))}
            unseen -= neighbors
            component |= neighbors
            todo.extend(neighbors)
        out.add(frozenset(component))
    return out


def independent_components(zs):
    """Union-find: for each coordinate, union its positive factorizations."""
    parent = list(range(len(zs)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for k in range(4):
        active = [i for i, z in enumerate(zs) if z[k] > 0]
        if active:
            for i in active[1:]:
                parent[root(i)] = root(active[0])
    groups = {}
    for i, z in enumerate(zs):
        groups.setdefault(root(i), set()).add(z)
    return {frozenset(c) for c in groups.values()}


def pure(i, n):
    return tuple(n if k == i else 0 for k in range(4))


def critical_coefficient(ds, i):
    others = ds[:i] + ds[i + 1:]
    # Each other generator supplies an explicit positive relation; no open search.
    upper = min(ds[j] // gcd(ds[i], ds[j]) for j in range(4) if j != i)
    return next(n for n in range(1, upper + 1) if factorizations(others, n * ds[i]))


def positive_coefficient(ds, m, i):
    assert m != i
    others = ds[:i] + ds[i + 1:]
    position = m if m < i else m - 1
    upper = ds[m] // gcd(ds[m], ds[i])
    return next(n for n in range(1, upper + 1)
                if any(z[position] > 0 for z in factorizations(others, n * ds[i])))


def is_minimal(ds):
    return all(not factorizations(ds[:i] + ds[i + 1:], ds[i]) for i in range(4))


def is_primary_from_critical_fibers(ds, a):
    # a_ii = b_(m,ii) iff the a_ii fiber has an i-free factorization positive at m.
    fibers = [factorizations(ds, a[i] * ds[i]) for i in range(4)]
    return any(sum(any(z[i] == 0 and z[m] > 0 for z in fibers[i])
                   for i in range(4) if i != m) >= 2 for m in range(4))


def fidelity():
    primary = (196, 225, 256, 289)
    a = tuple(critical_coefficient(primary, i) for i in range(4))
    assert a == (17, 7, 15, 4)
    assert is_minimal(primary) and gcd(*primary) == 1
    assert is_primary_from_critical_fibers(primary, a)
    secondary = (145, 203, 74, 111)
    a = tuple(critical_coefficient(secondary, i) for i in range(4))
    assert a == (7, 5, 3, 2)
    for i, pair in [(0, (0, 1)), (1, (0, 1)), (2, (2, 3)), (3, (2, 3))]:
        expected = {pure(j, a[j]) for j in pair}
        got = set(factorizations(secondary, a[i] * secondary[i]))
        assert got == expected == independent_factorizations(secondary, a[i] * secondary[i])
        assert len(graph_components(tuple(got))) == 2
    assert is_minimal(secondary) and gcd(*secondary) == 1
    assert not is_primary_from_critical_fibers(secondary, a)
    print('FIDELITY Example 2.9: primary a=(17,7,15,4); secondary a=(7,5,3,2);')
    print('  all four displayed secondary critical fibers reproduced exactly.')
    return secondary


def check_boundary(ds, i):
    p = gcd(ds[0], ds[1])
    b = positive_coefficient(ds, 0, i)
    r = p // gcd(p, ds[i])
    n = b * ds[i]
    zs = factorizations(ds, n)
    assert set(zs) == independent_factorizations(ds, n)
    components = graph_components(zs)
    assert components == independent_components(zs)
    assert b <= r
    assert (len(components) > 1) == (b == r), (ds, i, b, r, components)
    assert len(components) == (2 if b == r else 1), (ds, i, b, r, components)
    if b == r:
        left = frozenset(z for z in zs if z[2] == z[3] == 0)
        right = frozenset(z for z in zs if z[0] == z[1] == 0)
        assert left and right and components == {left, right}
    return b == r, b, r, len(zs), len(components)


def main():
    source_example = fidelity()
    for i in (2, 3):
        print(f'SOURCE boundary i={i}: {check_boundary(source_example, i)}')
    minimal = secondary = positive = negative = 0
    examples = {}
    pairings = [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]
    for ds in combinations(range(2, BOUND + 1), 4):
        if gcd(*ds) != 1 or not is_minimal(ds):
            continue
        minimal += 1
        for a0, a1, a2, a3 in pairings:
            left_degree = lcm(ds[a0], ds[a1])
            right_degree = lcm(ds[a2], ds[a3])
            if left_degree == right_degree:
                continue
            expected_left = {pure(a0, left_degree // ds[a0]),
                             pure(a1, left_degree // ds[a1])}
            expected_right = {pure(a2, right_degree // ds[a2]),
                              pure(a3, right_degree // ds[a3])}
            if set(factorizations(ds, left_degree)) != expected_left:
                continue
            if set(factorizations(ds, right_degree)) != expected_right:
                continue
            a = tuple(critical_coefficient(ds, i) for i in range(4))
            assert {a[a0] * ds[a0], a[a1] * ds[a1]} == {left_degree}
            assert {a[a2] * ds[a2], a[a3] * ds[a3]} == {right_degree}
            assert not is_primary_from_critical_fibers(ds, a)
            secondary += 1
            orders = [(a0, a1, a2, a3), (a1, a0, a2, a3),
                      (a0, a1, a3, a2), (a1, a0, a3, a2),
                      (a2, a3, a0, a1), (a3, a2, a0, a1),
                      (a2, a3, a1, a0), (a3, a2, a1, a0)]
            for order in orders:
                ordered = tuple(ds[k] for k in order)
                for i in (2, 3):
                    result = check_boundary(ordered, i)
                    examples.setdefault(result[0], (ordered, i, result))
                    if result[0]:
                        positive += 1
                    else:
                        negative += 1
    assert minimal == 28106 and secondary == 85
    assert positive == 196 and negative == 1164
    print(f'EXHAUSTIVE generators in [2,{BOUND}]: {minimal} minimally generated gcd-one unordered tuples;')
    print(f'  {secondary} secondary unordered tuples; {secondary * 8} source-compatible ordered tuples;')
    print(f'  {positive + negative} boundary cases: {positive} Betti, {negative} connected; no failures.')
    print('  Every checked fiber and component partition agrees with a second implementation.')
    for key in (True, False):
        print(f'EXAMPLE Betti={key}: {examples[key]}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
