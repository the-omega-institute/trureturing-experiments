#!/usr/bin/env python3
"""Exact tests of Doust--Weston section 5 and the proposed cut identity.

Every length is its rational coefficient of pi. No floating point, numerical
quadrature, or truncated infinite tail is used. These checks are finite tests,
not a proof of the universal countable-closed assertion.
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product


def circle_distance(a, b):
    d = (a - b) % 2
    return min(d, 2 - d)


def l1(a, b):
    return sum((abs(x - y) for x, y in zip(a, b)), Q(0))


def source_fidelity():
    # Source: x_k has longitude pi/2^k, its antipode adds pi,
    # and z is the north pole. gamma(x_k)_i = 0 for i<k,
    # 2^(-i-1) otherwise; antipodes complement the geometric series.
    max_k = 12
    points = [(k, antipodal) for k in range(max_k + 1)
              for antipodal in (False, True)] + ["north"]
    end = max_k + 1

    def coefficient(p, i):
        if p == "north":
            return Q(1, 2)
        k, antipodal = p
        base = Q(int(i >= k))
        return 1 - base if antipodal else base

    def exact_infinite_l1(p, q):
        prefix = sum((abs(coefficient(p, i) - coefficient(q, i))
                      * Q(1, 2 ** (i + 1)) for i in range(end)), Q(0))
        # From end onward both coefficient sequences are constant, and
        # sum_{i=end}^infinity 2^(-i-1) = 2^(-end) exactly.
        tail = abs(coefficient(p, end) - coefficient(q, end)) * Q(1, 2 ** end)
        return prefix + tail

    def expected(p, q):
        if p == q:
            return Q(0)
        if "north" in (p, q):
            return Q(1, 2)
        k, a = p
        j, b = q
        return circle_distance(Q(1, 2 ** k) + int(a),
                               Q(1, 2 ** j) + int(b))

    pairs = 0
    for p, q in product(points, repeat=2):
        assert exact_infinite_l1(p, q) == expected(p, q), (p, q)
        pairs += 1
    assert exact_infinite_l1((0, False), (0, True)) == 1
    print(f"FIDELITY: section 5 dyadic antipodal pairs plus north pole; "
          f"k=0..{max_k}; {pairs} ordered pairs; infinite tails summed exactly")


def positive_cos_pi(t):
    t %= 2
    assert t not in (Q(1, 2), Q(3, 2)), "representative is on an equator"
    return t < Q(1, 2) or t > Q(3, 2)


def hemisphere_cells(n):
    # K is the regular 2n-gon on the equator plus the two poles.
    # Its dual great circles have longitude boundaries 1/2+j/n (in pi units),
    # with the normal sphere's equator splitting each sector in two.
    # Normalized area of each half-lune is (pi/n)/(4*pi) = 1/(4n).
    points = [("equator", Q(j, n)) for j in range(2 * n)]
    points += [("pole", 1), ("pole", -1)]
    cells = []
    for j in range(2 * n):
        for upper in (False, True):
            def pattern(offset):
                longitude = Q(1, 2) + Q(j, n) + offset
                return tuple(positive_cos_pi(longitude - p[1]) if p[0] == "equator"
                             else upper == (p[1] == 1) for p in points)
            mid = pattern(Q(1, 2 * n))
            assert mid == pattern(Q(1, 4 * n)) == pattern(Q(3, 4 * n))
            cells.append((mid, Q(1, 4 * n)))
    assert len({pattern for pattern, _ in cells}) == 4 * n
    assert sum((w for _, w in cells), Q(0)) == 1
    return points, cells


def geodesic_pi(p, q):
    if p[0] == q[0] == "equator":
        return circle_distance(p[1], q[1])
    if p[0] != q[0]:
        return Q(1, 2)
    return Q(0 if p[1] == q[1] else 1)


def check_cut_embedding(points, cells, selected):
    # Aggregating equal restricted patterns models exactly the proof's
    # sign-pattern fibers, including disconnected fibers and the empty K.
    fibers = defaultdict(Q)
    for pattern, weight in cells:
        fibers[tuple(pattern[i] for i in selected)] += weight
    assert sum(fibers.values(), Q(0)) == 1
    vectors = [tuple(weight * int(pattern[j]) for pattern, weight in fibers.items())
               for j in range(len(selected))]
    for vector in vectors:
        assert sum(map(abs, vector), Q(0)) <= 1
    pairs = 0
    for i, j in product(range(len(selected)), repeat=2):
        expected = geodesic_pi(points[selected[i]], points[selected[j]])
        actual = l1(vectors[i], vectors[j])
        assert actual == expected, (points[selected[i]], points[selected[j]], actual, expected)
        assert (vectors[i] == vectors[j]) == (i == j)
        pairs += 1
    return pairs


def settlement_tests():
    pairs = 0
    for n in range(1, 17):
        points, cells = hemisphere_cells(n)
        pairs += check_cut_embedding(points, cells, tuple(range(len(points))))
    print(f"CUT_IDENTITY: regular 2n equatorial points plus both poles; "
          f"n=1..16; {pairs} ordered pairs; exact half-lune weights")
    points, cells = hemisphere_cells(2)
    subset_pairs = 0
    for mask in range(1 << len(points)):
        selected = tuple(i for i in range(len(points)) if mask & (1 << i))
        subset_pairs += check_cut_embedding(points, cells, selected)
    print(f"DEGENERACIES: all 64 subsets of the six-point octahedron, "
          f"including empty/singleton/antipodal sets; {subset_pairs} ordered pairs")
    # Negative control: a factor-of-two normalization defect must be detected.
    p, q = points.index(("pole", 1)), points.index(("pole", -1))
    wrong_distance = sum((weight / 2 for pattern, weight in cells if pattern[p] != pattern[q]), Q(0))
    assert wrong_distance != geodesic_pi(points[p], points[q])
    print("NEGATIVE_CONTROL: half-scale cut weights fail the antipodal distance")


if __name__ == "__main__":
    source_fidelity()
    settlement_tests()
    print("ALL_CHECKS_OK")
