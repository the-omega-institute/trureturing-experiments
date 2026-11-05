#!/usr/bin/env python3
"""Exact checks for arXiv:2610.19193v1 and the proposed Conjecture 2 settlement.

Standard-library Python 3 only. Coordinates, areas, sample points and polynomial
coefficients use integers or Fraction. Finite checks are not a proof of cofinality.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import isqrt


def seven(W, H, a, b, c, d, e, f):
    return [(a, W, 0, d), (0, a, 0, f), (0, c, f, H),
            (c, W, e, H), (b, W, d, e), (a, b, d, f), (b, c, e, f)]


def published(p, q):
    assert 0 < q < p < 7*q and p % 2 == q % 2
    W = H = 105*q
    cuts = (7*(7*q-p), 42*q, 7*(7*q+p), 5*(8*q-p),
            5*(13*q-p), Q(15, 2)*(7*q+p))
    return W, H, seven(W, H, *cuts)


def compressed(p, q):
    assert 0 < q < p < 7*q
    W, H = 15*q, 42*q
    cuts = (7*q-p, 6*q, 7*q+p, 2*(8*q-p), 2*(13*q-p), 3*(7*q+p))
    a, b, c, d, e, f = cuts
    assert 0 < a < b < c < W and 0 < d < e < f < H
    return W, H, seven(W, H, *cuts)


def check_partition(W, H, rects):
    assert len(rects) == 7
    dimensions = []
    for x0, x1, y0, y1 in rects:
        assert 0 <= x0 < x1 <= W and 0 <= y0 < y1 <= H
        assert all(Q(t).denominator == 1 for t in (x0, x1, y0, y1))
        dimensions.append((x1-x0, y1-y0))
    for r, s in combinations(rects, 2):
        assert max(r[0], s[0]) >= min(r[1], s[1]) or \
               max(r[2], s[2]) >= min(r[3], s[3])
    # Independent exact grid check of every open elementary cell, edge and vertex.
    # A rectangle's membership is constant in each elementary cell/edge.
    xs = sorted({0, W} | {r[i] for r in rects for i in (0, 1)})
    ys = sorted({0, H} | {r[i] for r in rects for i in (2, 3)})
    xprobe = xs + [Q(a+b, 2) for a, b in zip(xs, xs[1:])]
    yprobe = ys + [Q(a+b, 2) for a, b in zip(ys, ys[1:])]
    for x, y in product(xprobe, yprobe):
        assert any(x0 <= x <= x1 and y0 <= y <= y1 for x0, x1, y0, y1 in rects)
    shapes = [tuple(sorted(d)) for d in dimensions]
    assert len(set(shapes)) == 7  # includes rotations
    areas = [int(w*h) for w, h in dimensions]
    perimeters = [int(2*(w+h)) for w, h in dimensions]
    assert sum(areas) == W*H
    return dimensions, areas, perimeters, max(areas)-min(areas)


def next_pair(p, q):
    return 170*p+741*q, 39*p+170*q


def product_coefficients(x, y):
    """(a*p+b*q)(c*p+d*q) -> coefficients of p^2, p*q, q^2."""
    a, b = x
    c, d = y
    return a*c, a*d+b*c, b*d


def symbolic_checks():
    p2 = product_coefficients((170, 741), (170, 741))
    q2 = product_coefficients((39, 170), (39, 170))
    assert tuple(a-19*b for a, b in zip(p2, q2)) == (1, 0, -19)
    print('SYMBOLIC norm recurrence coefficients = (1, 0, -19)')
    dims = [((1, 8), (-2, 16)), ((-1, 7), (3, 21)),
            ((1, 7), (-3, 21)), ((-1, 8), (2, 16)),
            ((0, 9), (0, 10)), ((1, -1), (5, 5)), ((1, 1), (5, -5))]
    expected = [-2, -3, -3, -2, 0, 5, 5]
    perimeters = []
    for (w, h), off in zip(dims, expected):
        a, b, c = product_coefficients(w, h)
        assert b == 0
        assert (19*a+c, a) == (90, off)  # replace p^2 with 19*q^2+1
        perimeters.append(tuple(2*(a+b) for a, b in zip(w, h)))
    differences = [tuple(a-b for a, b in zip(perimeters[i], perimeters[j]))
                   for i, j in [(0, 3), (1, 2), (5, 6)]]
    assert differences == [(-4, 0), (8, 0), (0, 16)]
    print('SYMBOLIC seven area offsets from 90*q^2 =', expected)
    print('SYMBOLIC equal-area perimeter differences = -4*p, 8*p, 16*q')


def source_fidelity():
    W, H, rects = published(13, 3)
    dims, areas, perimeters, defect = check_partition(W, H, rects)
    expected = [(259, 55), (56, 255), (238, 60), (77, 185),
                (189, 75), (70, 200), (112, 125)]
    assert dims == expected and min(areas) == 14000 and max(areas) == 14280
    assert W == H == 315 and defect == 280
    print('FIDELITY source Example 1 canvas = 315 x 315')
    print('FIDELITY dimensions =', expected)
    print('FIDELITY min area = 14000; max area = 14280; defect = 280')
    p, q = 13, 3
    sides = []
    for j in range(8):
        W, H, rects = published(p, q)
        assert p*p-19*q*q == -2
        assert check_partition(W, H, rects)[3] == 280
        assert W % 630 == 315
        sides.append(W)
        p, q = next_pair(p, q)
    assert sides[:4] == [315, 106785, 36306585, 12344132115]
    assert all(sides[i+2] == 340*sides[i+1]-sides[i] for i in range(6))
    print('FIDELITY source square family checked at j = 0..7; table and recurrence agree')


def settlement_checks():
    p, q = 170, 39
    for j in range(32):
        assert p*p-19*q*q == 1 and q >= 39+j and q < p < 7*q
        W, H, rects = compressed(p, q)
        dims, areas, perimeters, defect = check_partition(W, H, rects)
        assert defect == 8 and H*5 == W*14
        C = 90*q*q
        assert areas == [C-2, C-3, C-3, C-2, C, C+5, C+5]
        assert [perimeters[0]-perimeters[3], perimeters[1]-perimeters[2],
                perimeters[5]-perimeters[6]] == [-4*p, 8*p, 16*q]
        SW, SH, source = published(2*p, 2*q)
        scaled = [(Q(x0, 14), Q(x1, 14), Q(y0, 5), Q(y1, 5))
                  for x0, x1, y0, y1 in source]
        assert (Q(SW, 14), Q(SH, 5), scaled) == (W, H, rects)
        assert check_partition(SW, SH, source)[3] == 560
        if j == 0:
            assert (W, H) == (585, 1638)
            assert areas == [136888, 136887, 136887, 136888, 136890, 136895, 136895]
            assert perimeters == [1532, 2864, 1504, 2212, 1482, 2352, 1728]
            print('SETTLEMENT first canvas =', W, 'x', H)
            print('SETTLEMENT first dimensions =', dims)
            print('SETTLEMENT first areas =', areas)
            print('SETTLEMENT first perimeters =', perimeters)
            print('SETTLEMENT first area sum =', sum(areas), '; defect =', defect)
        p1, q1 = next_pair(p, q)
        assert q1 >= q+1
        if j == 0:
            assert (p1, q1) == (57799, 13260)
        p, q = p1, q1
    print('SETTLEMENT 32 recurrence instances pass geometry, areas, congruence and affine checks')
    small = []
    for q in range(1, 1001):
        p = isqrt(19*q*q+1)
        if p*p == 19*q*q+1:
            assert check_partition(*compressed(p, q))[3] == 8
            small.append((p, q))
    assert small == [(170, 39)]
    print('SETTLEMENT exhaustive norm-1 solutions for 1 <= q <= 1000 =', small)


if __name__ == '__main__':
    symbolic_checks()
    source_fidelity()
    settlement_checks()
    print('ALL_CHECKS_OK')
