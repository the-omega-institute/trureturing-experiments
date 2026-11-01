#!/usr/bin/env python3
"""Exact finite checks for the LEDM branches of arXiv:2610.11950v1, Section 6.

This is finite evidence, not a proof of the unbounded claim. Rectangles are
labelled: combinations_with_replacement deliberately retains multiplicity.
Only Python integers, finite sets and bit masks are used.
"""
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product


def shape(a):
    return len(a), len(a[0]) if a else 0


def support(a):
    nr, nc = shape(a)
    assert all(a[i][j] >= 0 for i in range(nr) for j in range(nc))
    return frozenset((i, j) for i in range(nr) for j in range(nc) if a[i][j] > 0)


def strict_pairs(a):
    return tuple((p, q) for p, q in combinations(sorted(support(a)), 2)
                 if a[p[0]][p[1]] * a[q[0]][q[1]] > a[p[0]][q[1]] * a[q[0]][p[1]])


def valid(a, rects, refined):
    supp = support(a)
    for rect in rects:
        rows = {p[0] for p in rect}
        cols = {p[1] for p in rect}
        if rect != frozenset(product(rows, cols)) or not rect <= supp:
            return False
    if frozenset().union(*rects) != supp:
        return False
    return not refined or all(sum(p in rect or q in rect for rect in rects) >= 2
                              for p, q in strict_pairs(a))


def rectangles(a, maximal):
    nr, nc = shape(a)
    supp = support(a)
    found = set()
    for ib in range(1, 1 << nr):
        rows = tuple(i for i in range(nr) if ib >> i & 1)
        if maximal:
            cols = tuple(j for j in range(nc) if all((i, j) in supp for i in rows))
            if cols:
                closed_rows = tuple(i for i in range(nr) if all((i, j) in supp for j in cols))
                found.add(frozenset(product(closed_rows, cols)))
        else:
            for jb in range(1, 1 << nc):
                cols = tuple(j for j in range(nc) if jb >> j & 1)
                rect = frozenset(product(rows, cols))
                if rect <= supp:
                    found.add(rect)
    return tuple(sorted(found, key=lambda s: tuple(sorted(s))))


@lru_cache(None)
def minimum(a, refined=True, maximal=True):
    cells = tuple(sorted(support(a)))
    if not cells:
        return 0, ()
    pairs = strict_pairs(a) if refined else ()
    rects = rectangles(a, maximal)
    effects = []
    for rect in rects:
        flags = [p in rect for p in cells]
        flags += [p in rect or q in rect for p, q in pairs]
        effects.append(sum(int(v) << i for i, v in enumerate(flags)))
    need_once = (1 << len(cells)) - 1
    need_twice = ((1 << len(pairs)) - 1) << len(cells)
    # Enlarging each rectangle preserves coverage and cannot remove a touching
    # label, so restricting the minimization to maximal rectangles is exact.
    for k in range(1, len(cells) + 1):
        for labels in combinations_with_replacement(range(len(rects)), k):
            once = twice = 0
            for label in labels:
                effect = effects[label]
                twice |= once & effect
                once |= effect
            if once & need_once == need_once and twice & need_twice == need_twice:
                witness = tuple(rects[i] for i in labels)
                assert valid(a, witness, refined)  # independent set-based check
                return k, witness
    raise AssertionError("singleton upper bound failed")


def restrict(a, rows, cols):
    return tuple(tuple(a[i][j] for j in cols) for i in rows)


def subsets(n):
    return [tuple(i for i in range(n) if bits >> i & 1) for bits in range(1 << n)]


def check_restrictions(a, rects, refined=True):
    assert valid(a, rects, refined)
    nr, nc = shape(a)
    count = 0
    for rows in subsets(nr):
        for cols in subsets(nc):
            b = restrict(a, rows, cols)
            rr = tuple(frozenset((ii, jj) for ii, i in enumerate(rows)
                                for jj, j in enumerate(cols) if (i, j) in rect)
                       for rect in rects)
            assert valid(b, rr, refined)
            if refined:
                for p, q in strict_pairs(b):
                    oldp = rows[p[0]], cols[p[1]]
                    oldq = rows[q[0]], cols[q[1]]
                    old = {l for l, rect in enumerate(rects) if oldp in rect or oldq in rect}
                    new = {l for l, rect in enumerate(rr) if p in rect or q in rect}
                    assert old == new
            count += 1
    return count


def ledm(coords):
    return tuple(tuple((i - j) ** 2 for j in coords) for i in coords)


def main():
    example1 = ((1, 1, 1), (1, 2, 2), (1, 2, 3))
    example2 = ((1, 2, 0), (1, 1, 1), (0, 1, 1))
    for label, a, expected in [("Example 1", example1, (1, 2)),
                                ("Example 2", example2, (2, 3))]:
        observed = minimum(a, False)[0], minimum(a, True)[0]
        assert observed == expected
        print(f"FIDELITY {label}: RCB={observed[0]} RRCB={observed[1]}")
    assert minimum(example1, True)[1][0] == minimum(example1, True)[1][1]
    print("FIDELITY Example 1: two coincident rectangles keep two distinct labels")

    rcb_values, rrcb_values = [], []
    restrictions = 0
    for n in range(1, 6):
        a = ledm(range(1, n + 1))
        k0, ordinary = minimum(a, False)
        k1, refined = minimum(a, True)
        rcb_values.append(k0)
        rrcb_values.append(k1)
        restrictions += check_restrictions(a, ordinary, False)
        restrictions += check_restrictions(a, refined, True)
        print(f"SMALL LEDM n={n}: RCB={k0} RRCB={k1}")
    for n in range(1, 6):
        for m in range(n, 6):
            assert rcb_values[n-1] <= rcb_values[m-1]
            assert rrcb_values[n-1] <= rrcb_values[m-1]
            assert ledm(range(1, n+1)) == restrict(ledm(range(1, m+1)), tuple(range(n)), tuple(range(n)))
    print("SMALL monotonicity: every 1 <= n <= m <= 5, both LEDM branches")

    raw_families = 0
    for entries in product(range(3), repeat=4):
        a = (entries[:2], entries[2:])
        for refined in [False, True]:
            assert minimum(a, refined, True)[0] == minimum(a, refined, False)[0]
            restrictions += check_restrictions(a, minimum(a, refined)[1], refined)
        raw_rects = (frozenset(),) + rectangles(a, False)
        for k in range(3):
            for fam in combinations_with_replacement(raw_rects, k):
                if valid(a, fam, True):
                    restrictions += check_restrictions(a, fam, True)
                    raw_families += 1
    print(f"EXHAUSTIVE all 81 nonnegative 2x2 matrices with entries 0,1,2: maximal/raw minima agree")
    print(f"EXHAUSTIVE {raw_families} refined families of size <= 2, including empty/repeated rectangles")

    a = ledm((-2, -2, 0, 3, 3))
    restrictions += check_restrictions(a, minimum(a, True)[1], True)
    print("BOUNDARIES zero matrices, empty restrictions, one-cell matrices, repeated real coordinates checked")
    print(f"RESTRICTIONS {restrictions} labelled family/submatrix checks; exact touching-label sets preserved")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
