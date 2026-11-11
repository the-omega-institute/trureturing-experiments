#!/usr/bin/env python3
"""Exact finite checks for Teheux, arXiv:2611.12469v1, Question 8.4(2).

Only integers, finite sets and bit masks are used. No external packages.
These finite checks do not establish the assertion on infinite structures.
"""
from itertools import combinations, combinations_with_replacement, product


def majority(a, b, c):
    return (a & b) | (a & c) | (b & c)


def star_configurations(ultrafilters):
    """Literal sigma: five distinct points, all four leaf triples centred at U."""
    for centre in ultrafilters:
        for leaves in combinations([v for v in ultrafilters if v != centre], 4):
            if all(majority(*triple) == centre for triple in combinations(leaves, 3)):
                yield centre, leaves


def tree_fidelity(n, edges):
    neighbours = [set() for _ in range(n)]
    for a, b in edges:
        neighbours[a].add(b)
        neighbours[b].add(a)

    def component(a, b):
        seen = {a}
        todo = [a]
        while todo:
            x = todo.pop()
            for y in neighbours[x]:
                if {x, y} == {a, b} or y in seen:
                    continue
                seen.add(y)
                todo.append(y)
        return sum(1 << x for x in seen)

    whole = (1 << n) - 1
    halves = [0, whole]
    for a, b in edges:
        h = component(a, b)
        halves += [h, whole ^ h]
    assert len(set(halves)) == 2 + 2 * len(edges)
    up = [sum(1 << j for j, t in enumerate(halves) if s & t == s) for s in halves]
    complement = [halves.index(whole ^ s) for s in halves]
    ufs = ultrafilters_of(up, complement)
    principal = [sum(1 << i for i, s in enumerate(halves) if s >> v & 1) for v in range(n)]
    assert set(ufs) == set(principal)

    def path(a, b):
        todo = [(a, [a])]
        while todo:
            x, route = todo.pop()
            if x == b:
                return set(route)
            todo += [(y, route + [y]) for y in neighbours[x] if y not in route]
        raise AssertionError('disconnected source tree')

    for a, b, c in product(range(n), repeat=3):
        intersection = path(a, b) & path(a, c) & path(b, c)
        assert len(intersection) == 1
        v = next(iter(intersection))
        assert majority(principal[a], principal[b], principal[c]) == principal[v]
    sigma = not any(star_configurations(ufs))
    assert sigma == first_order(up, complement)
    return sigma, len(halves), len(ufs)


def ultrafilters_of(up, complement):
    n = len(up)
    out = []
    for mask in range(1, (1 << n) - 1):
        if not all(bool(mask >> a & 1) != bool(mask >> complement[a] & 1) for a in range(n)):
            continue
        if all(not (mask >> a & 1) or up[a] & mask == up[a] for a in range(n)):
            out.append(mask)
    return out


def obstruction(up, complement, ps):
    c = 0
    for p in ps:
        c |= up[complement[p]]
    for p, q in combinations(ps, 2):
        c |= up[p] & up[q]
    return c


def consistent(mask, complement):
    return not any((mask >> a & 1) and (mask >> complement[a] & 1) for a in range(len(complement)))


def first_order(up, complement):
    """Literal forall p1...p4: a zero parameter or exists x in C with x* in C.

    Zero parameters make the disjunction true. Each nonzero ordered tuple
    is represented by its nondecreasing rearrangement, since C is symmetric.
    Repetitions are included, with the six indexed pair clauses retained.
    """
    return all(not consistent(obstruction(up, complement, ps), complement)
               for ps in combinations_with_replacement(range(1, len(up)), 4))


def poc_sets(r):
    """Every bounded poc set on 2+2r points, up to complement-preserving labels.

    Bottom=0, top=1; complement is a xor 1. For two proper complement
    pairs there are four possible dual orbits of cross-comparisons. At most
    one is allowed: two force either a nonzero a<=a* or a cycle. Enumerate
    the five choices (none or one orbit), then retain exactly transitive
    relations. No transitive closure is taken, so each relation occurs once.
    """
    n = 2 + 2 * r
    orbits = []
    for a, b in combinations(range(2, n, 2), 2):
        orbits.append([(a, b), (a, b ^ 1), (a ^ 1, b), (a ^ 1, b ^ 1)])
    for choices in product(range(5), repeat=len(orbits)):
        up = [1 << a | 2 for a in range(n)]
        up[0] = (1 << n) - 1
        for options, choice in zip(orbits, choices):
            if choice:
                a, b = options[choice - 1]
                up[a] |= 1 << b
                up[b ^ 1] |= 1 << (a ^ 1)
        if any(up[b] & up[a] != up[b] for a in range(n) for b in range(n) if up[a] >> b & 1):
            continue
        complement = [a ^ 1 for a in range(n)]
        # Independently validate every source axiom on the retained relation.
        assert all(up[a] >> a & 1 for a in range(n))
        assert all(a == b or not (up[a] >> b & 1 and up[b] >> a & 1)
                   for a, b in product(range(n), repeat=2))
        assert all(bool(up[a] >> b & 1) == bool(up[b ^ 1] >> (a ^ 1) & 1)
                   for a, b in product(range(n), repeat=2))
        assert all(bool(up[a] >> (a ^ 1) & 1) == (a == 0) for a in range(n))
        yield up, complement


def check_frame(up, complement):
    n = len(up)
    ufs = ultrafilters_of(up, complement)
    assert ufs
    assert all(majority(a, b, c) in ufs for a, b, c in product(ufs, repeat=3))
    sigma = not any(star_configurations(ufs))
    assert sigma == first_order(up, complement), (up, complement, sigma)
    # Direct ordered quantifiers give an independent check of the symmetry reduction.
    if n <= 8:
        literal = all(not consistent(obstruction(up, complement, ps), complement)
                      for ps in product(range(1, n), repeat=4))
        assert literal == sigma
    extensions = 0
    for mask in range(1, 1 << n):
        if consistent(mask, complement) and all(not (mask >> a & 1) or up[a] & mask == up[a]
                                               for a in range(n)):
            assert any(mask & u == mask for u in ufs)
            extensions += 1
    constructions = 0
    for ps in combinations_with_replacement(range(1, n), 4):
        c = obstruction(up, complement, ps)
        if not consistent(c, complement):
            continue
        assert c and len(set(ps)) == 4
        centres = [u for u in ufs if c & u == c]
        assert centres
        for u in centres:
            leaves = []
            for p in ps:
                down = sum(1 << a for a in range(n) if up[a] >> complement[p] & 1)
                v = (u & ~down) | up[p]
                assert v in ufs and v != u
                leaves.append(v)
            assert len(set([u] + leaves)) == 5
            assert all(majority(*t) == u for t in combinations(leaves, 3))
            assert all(sum((v ^ u) >> a & 1 for v in leaves) <= 1 for a in range(n))
            constructions += 1
    # Verify the forward extraction, independently enumerating all stars.
    stars = 0
    for u, leaves in star_configurations(ufs):
        witnesses = [next(p for p in range(n) if v >> p & 1 and not (u >> p & 1)) for v in leaves]
        c = obstruction(up, complement, witnesses)
        assert all(witnesses) and c & u == c and consistent(c, complement)
        stars += 1
    return sigma, extensions, constructions, stars


def main():
    d = tree_fidelity(6, [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)])
    s = tree_fidelity(5, [(0, 1), (0, 2), (0, 3), (0, 4)])
    assert d == (True, 12, 6) and s == (False, 10, 5)
    print('FIDELITY Remark 7.5: double trivalent tree sigma=True; four-leaf-star quotient sigma=False')
    print('FIDELITY tree median = majority on all vertex triples; principal ultrafilters exhaust both complex algebras')
    totals = [0, 0, 0, 0, 0, 0]
    for r in range(5):
        count = good = extensions = constructions = stars = 0
        for up, complement in poc_sets(r):
            sigma, e, c, t = check_frame(up, complement)
            count += 1
            good += sigma
            extensions += e
            constructions += c
            stars += t
        row = [count, good, count - good, extensions, constructions, stars]
        totals = [a + b for a, b in zip(totals, row)]
        print(f'P_SIZE={2+2*r} frames={count} sigma_true={good} sigma_false={count-good} '
              f'consistent_upsets={extensions} flip_constructions={constructions} extracted_stars={stars}', flush=True)
    print('TOTAL frames=%d sigma_true=%d sigma_false=%d consistent_upsets=%d flip_constructions=%d extracted_stars=%d' % tuple(totals))
    print('COVERAGE all bounded poc sets of cardinality at most 10 up to complement-preserving relabelling; no infinite-case claim')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
