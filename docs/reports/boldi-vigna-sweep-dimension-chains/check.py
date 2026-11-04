#!/usr/bin/env python3
"""Exact finite checks for Boldi--Vigna, arXiv:2610.41385v1, Section 6.

Only Python's standard library is used. No floating point, Lean claim proof,
or numerical optimization is involved. Factor sizes, rather than top ranks,
are printed. Convexity and all separation obligations are exhaustively tested.
"""
from itertools import product, permutations
from math import comb


def grid(sizes):
    return list(product(*(range(n) for n in sizes)))


def le(x, y):
    return all(a <= b for a, b in zip(x, y))


def symmetric_chains(sizes):
    chains = [[()]]
    for size in sizes:
        b = size - 1
        extended = []
        for chain in chains:
            a = len(chain) - 1
            for q in range(min(a, b) + 1):
                extended.append([chain[q] + (j,) for j in range(b-q+1)]
                                + [chain[i] + (b-q,) for i in range(q+1, a+1)])
        chains = extended
    return chains


def rank_counts(points):
    R = max(map(sum, points))
    return [sum(sum(x) == j for x in points) for j in range(R+1)]


def convolution_counts(sizes):
    a = [1]
    for n in sizes:
        out = [0] * (len(a)+n-1)
        for i, ai in enumerate(a):
            for j in range(n):
                out[i+j] += ai
        a = out
    return a


def target(points):
    a = rank_counts(points)
    R = len(a)-1
    return 1 if R == 0 else a[R//2] + a[(R-1)//2]


def intervals(points, relation):
    n = len(points)
    between = [[sum(1 << z for z in range(n)
                    if relation(points[x], points[z]) and relation(points[z], points[y]))
                for y in range(n)] for x in range(n)]
    answer = []
    for mask in range(1 << n):
        members = [i for i in range(n) if mask >> i & 1]
        if all(between[x][y] & mask == between[x][y]
               for x in members for y in members):
            answer.append(mask)
    return answer, between


def obligations(convex, n):
    return [(mask, y) for mask in convex for y in range(n) if not mask >> y & 1]


def separation_mask(order, requirements, n):
    positions = [0]*n
    for j, x in enumerate(order):
        positions[x] = j
    mask = 0
    for i, (I, y) in enumerate(requirements):
        if all(positions[y] < positions[x] for x in range(n) if I >> x & 1):
            mask |= 1 << i
    return mask


def check_orders(orders, requirements, n):
    covered = 0
    for order in orders:
        assert sorted(order) == list(range(n))
        covered |= separation_mask(order, requirements, n)
    assert covered == (1 << len(requirements)) - 1


def exact_dimension_small(points, relation):
    """Exhaust all orders; solve the finite set-cover problem exactly (n <= 6)."""
    n = len(points)
    assert n <= 6
    convex, _ = intervals(points, relation)
    req = obligations(convex, n)
    cover = set(separation_mask(p, req, n) for p in permutations(range(n)))
    full = (1 << len(req))-1
    # Remove redundant equal cover masks; search by an uncovered obligation.
    supports = [[m for m in cover if m >> i & 1] for i in range(len(req))]
    def feasible(k):
        failed = set()
        def visit(left, remaining):
            if not left:
                return True
            if remaining == 0 or (left, remaining) in failed:
                return False
            bits = [i for i in range(len(req)) if left >> i & 1]
            i = min(bits, key=lambda j: len(supports[j]))
            options = sorted({m & left for m in supports[i]}, key=int.bit_count, reverse=True)
            for m in options:
                if visit(left & ~m, remaining-1):
                    return True
            failed.add((left, remaining))
            return False
        return visit(full, k)
    for k in range(n+1):
        if feasible(k):
            return k
    raise AssertionError('singleton-first orders must suffice')


def check_grid(sizes, exhaustive=True):
    points = grid(sizes)
    chains = symmetric_chains(sizes)
    R = sum(n-1 for n in sizes)
    flat = [x for c in chains for x in c]
    assert len(flat) == len(set(flat)) == len(points)
    assert set(flat) == set(points)
    for c in chains:
        assert sum(c[0]) + sum(c[-1]) == R
        assert all(le(x, y) and sum(y)-sum(x) == 1 for x, y in zip(c, c[1:]))
    a = rank_counts(points)
    assert a == convolution_counts(sizes) == a[::-1]
    cost = sum(min(2, len(c)) for c in chains)
    assert cost == target(points)
    if not exhaustive:
        return cost
    indices = {x: i for i, x in enumerate(points)}
    orders = []
    for c in chains:
        rest = [x for x in points if x not in c]
        orders.append([indices[x] for x in c+rest])
        if len(c) > 1:
            orders.append([indices[x] for x in c[::-1]+rest])
    convex, between = intervals(points, le)
    check_orders(orders, obligations(convex, len(points)), len(points))
    if R == 0:
        # Empty interval and the sole outsider require a nonempty tuple.
        assert (0, 0) in obligations(convex, 1) and cost == 1
    else:
        ranks = {R//2, (R-1)//2} if R % 2 == 0 else {R//2, R//2+1}
        S = [i for i, x in enumerate(points) if sum(x) in ranks]
        assert len(S) == cost
        allS = sum(1 << i for i in S)
        witnesses = {}
        for y in S:
            others = [i for i in S if i != y]
            hull = 0
            for u in others:
                for v in others:
                    hull |= between[u][v]
            assert hull in convex and not hull >> y & 1
            assert hull & allS == allS & ~(1 << y)
            witnesses[y] = hull
        # A sweep serving both y and z would require y < z and z < y.
        assert all(witnesses[y] >> z & 1 and witnesses[z] >> y & 1
                   for y in S for z in S if y != z)
    if len(points) <= 6:
        assert exact_dimension_small(points, le) == cost
    print(f'GRID sizes={sizes} points={len(points)} convex={len(convex)} sweepdim={cost}')
    return cost


def main():
    # Source fidelity: Proposition 3 and the computations stated after it.
    chain_values = []
    for n in range(1, 7):
        d = exact_dimension_small(list(range(n)), lambda x, y: x <= y)
        assert d == (1 if n == 1 else 2)
        chain_values.append(d)
    print(f'FIDELITY chains n=1..6 dimensions={chain_values}')
    anti_values = []
    for n in range(1, 6):
        d = exact_dimension_small(list(range(n)), lambda x, y: x == y)
        assert d == n
        anti_values.append(d)
    print(f'FIDELITY antichains n=1..5 dimensions={anti_values}')
    for lower, upper in [(1, 1), (1, 2), (2, 2), (2, 3)]:
        pts = list(range(lower+upper))
        rel = lambda x, y: x == y or x < lower <= y
        d = exact_dimension_small(pts, rel)
        assert d == lower+upper
        print(f'FIDELITY complete_bipartite=({lower},{upper}) sweepdim={d}')
    examples = [(1,), (2,), (5,), (2, 2), (2, 3), (3, 3),
                (2, 2, 2), (2, 2, 3), (2, 2, 2, 2)]
    for sizes in examples:
        check_grid(sizes)
    # All ordered factor-size tuples in this stated finite range.
    cases = 0
    for d in range(1, 5):
        for sizes in product(range(1, 5), repeat=d):
            if len(grid(sizes)) <= 12 and sizes not in examples:
                check_grid(sizes)
                cases += 1
    print(f'EXHAUSTIVE_RANGE dimensions=1..4 factor_sizes=1..4 points<=12 additional_cases={cases}')
    for n in range(13):
        sizes = (2,)*n
        points = grid(sizes)
        a = rank_counts(points)
        assert a == [comb(n, j) for j in range(n+1)]
        value = check_grid(sizes, exhaustive=False)
        assert value == comb(n+1, (n+1)//2)
        print(f'BOOLEAN n={n} points={len(points)} sweep_count={value}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
