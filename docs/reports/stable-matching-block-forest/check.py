#!/usr/bin/env python3
"""Exact finite checks for Ishida, arXiv:2610.05617v1, Problem D.

Independent implementation from the mathematical definitions and the proposed
cyclic-cell construction. No upstream program is imported or copied.
Only integers, tuples, sets and exhaustive finite loops are used.
"""
from collections import Counter
from functools import cache
from itertools import permutations, product
from math import factorial


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def members(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def prefixes(rows):
    result = []
    for row in rows:
        seq = [0]
        for x in row:
            seq.append(seq[-1] | (1 << x))
        result.append(seq)
    return result


def blocks(P, R):
    """All paired blocks, allowing unequal index sets on the two sides.

    For each nonempty A, any block's B must be the top |A| set of its
    least-indexed proposer. Every member on both sides is then checked.
    """
    pp, rr = prefixes(P), prefixes(R)
    result = set()
    for A in range(1, 1 << len(P)):
        k = A.bit_count()
        p = (A & -A).bit_length() - 1
        B = pp[p][k]
        if all(pp[x][k] == B for x in members(A)) and all(
            rr[y][k] == A for y in members(B)):
            result.add((A, B))
    return result


def brute_blocks(P, R):
    """Separate rank-inequality implementation of the literal source definition."""
    n = len(P)
    pr = [{x: i for i, x in enumerate(row)} for row in P]
    rr = [{x: i for i, x in enumerate(row)} for row in R]
    result = set()
    for A in range(1, 1 << n):
        for B in range(1, 1 << n):
            if A.bit_count() != B.bit_count():
                continue
            if all(pr[p][r] < pr[p][s]
                   for p in members(A) for r in members(B)
                   for s in range(n) if not B & (1 << s)) and all(
                       rr[r][p] < rr[r][q]
                       for r in members(B) for p in members(A)
                       for q in range(n) if not A & (1 << q)):
                result.add((A, B))
    return result


def finest_shape(n, family):
    """Enumerate exact covers of both sides; select a cover with most parts."""
    full = (1 << n) - 1

    def visit(used_p, used_r, sizes):
        if used_p == full:
            if used_r == full:
                yield tuple(sorted(sizes))
            return
        first = (full ^ used_p) & -(full ^ used_p)
        for A, B in family:
            if A & first and not A & used_p and not B & used_r:
                yield from visit(used_p | A, used_r | B, sizes + (A.bit_count(),))

    candidates = tuple(visit(0, 0, ()))
    require(bool(candidates), 'the full block must be an exact cover')
    return max(candidates, key=lambda x: (len(x), x))


@cache
def shapes(n):
    """Every unordered size-labelled rooted tree, with child sum at most n.

    Child shapes are a multiset of previously generated shapes of size < n.
    Each multiset is generated once by nondecreasing indices into that pool.
    Uncovered agents are implicit, not additional tree nodes.
    """
    pool = tuple(s for k in range(1, n) for s in shapes(k))
    result = []

    def visit(start, remaining, children):
        result.append((n, children))
        for i in range(start, len(pool)):
            child = pool[i]
            if child[0] <= remaining:
                visit(i, remaining - child[0], children + (child,))

    visit(0, n, ())
    require(len(result) == len(set(result)), 'duplicate shape')
    return tuple(result)


def construct(shape, start=0, reverse_cells=False, reverse_labels=False):
    n, children = shape
    S = tuple(range(start, start + n))
    family = {frozenset(S)}
    if n == 1:
        require(not children, 'singleton has no proper child')
        return {start: (start,)}, {start: (start,)}, family
    cells = []
    P, R = {}, {}
    position = start
    for child in children:
        cp, cr, cf = construct(child, position, reverse_cells, reverse_labels)
        C = tuple(range(position, position + child[0]))
        cells.append((C, True))
        P.update(cp)
        R.update(cr)
        family.update(cf)
        position += child[0]
    for x in range(position, start + n):
        cells.append(((x,), False))
    if reverse_cells:
        cells.reverse()
    if reverse_labels:
        cells = [(tuple(reversed(C)), prescribed) for C, prescribed in cells]
    require(len(cells) >= 2, 'proper children force at least two cells')
    for j, (C, prescribed) in enumerate(cells):
        tail = tuple(x for offset in range(1, len(cells))
                     for x in cells[(j + offset) % len(cells)][0])
        for x in C:
            P[x] = P[x] + tail if prescribed else tail + (x,)
            R[x] = R[x] + tail if prescribed else (x,) + tail
    return P, R, family


def check_construction(shape, reverse_cells=False, reverse_labels=False):
    pd, rd, family = construct(shape, reverse_cells=reverse_cells,
                               reverse_labels=reverse_labels)
    n = shape[0]
    P = tuple(pd[i] for i in range(n))
    R = tuple(rd[i] for i in range(n))
    for row in P + R:
        require(len(row) == n and set(row) == set(range(n)), 'strict completeness')
    require(all(R[i][0] == i for i in range(n)), 'receiver diagonal first choices')
    expected = {(sum(1 << x for x in A),) * 2 for A in family}
    actual = blocks(P, R)
    require(actual == expected, ('extra/missing block', shape, expected, actual))
    if n <= 4:
        require(actual == brute_blocks(P, R), 'independent block detector disagreement')
    return P, R


def fidelity():
    """Theorem 9 / Table 3 and the stated size-three finest-shape census."""
    for n, expected in [(1, {(1,): 1}), (2, {(2,): 14, (1, 1): 2}),
                        (3, {(3,): 45768, (1, 2): 504, (1, 1, 1): 384})]:
        orders = tuple(permutations(range(n)))
        profiles = tuple(product(orders, repeat=n))
        counts = Counter()
        for P in profiles:
            for R in profiles:
                family = blocks(P, R)
                require(family == brute_blocks(P, R), 'literal prefix/rank disagreement')
                counts[finest_shape(n, family)] += 1
        require(dict(counts) == expected, ('Theorem 9 census', n, counts))
        require(sum(counts.values()) == factorial(n) ** (2 * n), 'profile total')
        print('THEOREM_9', 'n=' + str(n), 'profiles=' + str(sum(counts.values())),
              'finest_shapes=' + str(sorted(counts.items())))


def main():
    fidelity()
    source_counts = (1, 3, 10, 40, 170, 785, 3770)
    checked = 0
    for n, expected in enumerate(source_counts, 1):
        cases = shapes(n)
        require(len(cases) == expected, ('Problem D shape census', n, len(cases)))
        prime_compatible = sum(sum(c[0] for c in shape[1]) < n for shape in cases)
        if n >= 2:
            require(prime_compatible == 2 * len(shapes(n - 1)), 'prime-compatible census')
        for shape in cases:
            for reverse_cells, reverse_labels in product((False, True), repeat=2):
                check_construction(shape, reverse_cells, reverse_labels)
                checked += 1
        print('PROBLEM_D', 'n=' + str(n), 'shapes=' + str(len(cases)),
              'prime_compatible=' + str(prime_compatible), 'cell_order_variants=4',
              'all_static_blocks_exact=True')
    example = (5, ((2, ((1, ()),)), (2, ())))
    P, R = check_construction(example)
    require(P == ((0, 1, 2, 3, 4), (0, 1, 2, 3, 4), (3, 2, 4, 0, 1),
                  (2, 3, 4, 0, 1), (0, 1, 2, 3, 4)), 'supplied proposer example')
    require(R == ((0, 1, 2, 3, 4), (1, 0, 2, 3, 4), (2, 3, 4, 0, 1),
                  (3, 2, 4, 0, 1), (4, 0, 1, 2, 3)), 'supplied receiver example')
    print('SUPPLIED_EXAMPLE n=5 exact_blocks=' + str(sorted(blocks(P, R))))
    print('CONSTRUCTION_PROFILES_CHECKED=' + str(checked))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
