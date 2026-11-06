#!/usr/bin/env python3
"""Exact checks for arXiv:2610.14546v1 and an all-row replacement route.

The small models contain zero/negative entries: they test identities, not the
positive order-ten existence claim. No finite test establishes unboundedness.
Only Python's standard library is required. No source ancillary code is used.
"""
from fractions import Fraction
from itertools import combinations, product
from math import gcd, lcm


def det(a):
    """Bareiss elimination with row pivoting and checked exact divisions."""
    n = len(a)
    assert all(len(row) == n for row in a)
    if not n:
        return 1
    b = [list(row) for row in a]
    sign, previous = 1, 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if b[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            b[k], b[pivot] = b[pivot], b[k]
            sign = -sign
        p = b[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = b[i][j] * p - b[i][k] * b[k][j]
                assert numerator % previous == 0
                b[i][j] = numerator // previous
            b[i][k] = 0
        previous = p
    return sign * b[-1][-1]


def det_fraction(a):
    """Independent rational Gaussian elimination."""
    b = [[Fraction(v) for v in row] for row in a]
    answer = Fraction(1)
    for k in range(len(b)):
        pivot = next((i for i in range(k, len(b)) if b[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            b[k], b[pivot] = b[pivot], b[k]
            answer = -answer
        p = b[k][k]
        answer *= p
        for i in range(k + 1, len(b)):
            factor = b[i][k] / p
            for j in range(k, len(b)):
                b[i][j] -= factor * b[k][j]
    assert answer.denominator == 1
    return answer.numerator


def cube(a):
    return [[v**3 for v in row] for row in a]


def minor(a, i, j):
    return [[v for c, v in enumerate(row) if c != j]
            for r, row in enumerate(a) if r != i]


def cofactor(a, i, j):
    return (-1)**(i + j) * det(minor(a, i, j))


def gamma(e):
    return [(-1)**j * det([row[:j] + row[j + 1:] for row in e])
            for j in range(len(e) + 1)]


def rowcomb(u, e):
    return [sum(u[i] * e[i][j] for i in range(len(e)))
            for j in range(len(e[0]))]


def border(e, r, u, v, w, z):
    assert r == rowcomb(u, e)
    assert [a**3 for a in r] == rowcomb(v, cube(e))
    ep = [row + [wi] for row, wi in zip(e, w)] + [r + [z]]
    k1 = z - sum(ui * wi for ui, wi in zip(u, w))
    k3 = z**3 - sum(vi * wi**3 for vi, wi in zip(v, w))
    assert gamma(ep) == [k1 * a for a in gamma(e)] + [0]
    assert gamma(cube(ep)) == [k3 * a for a in gamma(cube(e))] + [0]
    return ep, k1, k3


def source_fidelity():
    # Source equation (2.1), entered directly from main.tex.
    b = [
        [4,3,5,5,4,3,5,4,2,3],
        [4,5,5,2,2,3,4,5,3,5],
        [5,5,3,4,5,4,5,3,5,3],
        [2,3,4,5,2,4,4,3,3,3],
        [4,3,3,2,3,4,2,2,5,2],
        [5,4,5,5,2,5,5,4,3,5],
        [2,4,2,3,5,3,5,4,5,2],
        [4,5,2,4,3,4,2,2,5,4],
        [4,5,2,2,3,4,3,5,4,4],
    ]
    expected_c = [4725,18,-1432,4682,-7509,-7230,504,4797,5531,-5416]
    expected_d = [340739472931722111,28733076063275820,
                  -369435421691654920,124984744924291244,
                  -1235085791262920481,-722001405221789442,
                  791954706125920512,183749257820706315,
                  484371424027128179,-300611788966121134]
    c, d = gamma(b), gamma(cube(b))
    assert c == expected_c and d == expected_d
    for e in [b, cube(b)]:
        for j in range(10):
            m = [row[:j] + row[j+1:] for row in e]
            assert det(m) == det_fraction(m)
    assert gcd(*c) == gcd(*d) == 1
    assert -159*c[0] + 41658*c[1] - c[2] == 1
    assert lcm(*(gcd(d[i], d[j]) for i, j in combinations(range(10), 2))) == 840
    r = b[0]
    D12 = c[0]*d[1]*r[1]**2 - c[1]*d[0]*r[0]**2
    assert D12 == 1123741091386468277532
    # Every duplicate-row point is on the homogeneous cone. It also witnesses
    # nonvanishing of EVERY lower-row cofactor polynomial on that cone.
    for i, r in enumerate(b):
        a = [r] + b
        assert det(a) == det(cube(a)) == 0
        for j in range(10):
            assert cofactor(a, i+1, j) == -c[j] != 0
            assert cofactor(cube(a), i+1, j) == -d[j] != 0
    print('FIDELITY: all 20 seed cofactors agree with two exact algorithms')
    print('FIDELITY: equations (3.3)-(3.5) and (3.8) reproduced; D12 =', D12)
    print('CONE: all 180 lower-cofactor nonvanishing witnesses verified')


def all_small_replacement_identities():
    # All 3^6 two-by-three lower blocks and all 5^3 possible top rows.
    # This verifies degree-1/degree-3 scaling for lower-row cofactors, Laplace
    # identities, and the precise singularity criterion whenever both dets=1.
    count, simultaneous, replacements = 0, 0, 0
    for entries in product(range(-1, 2), repeat=6):
        e = [list(entries[:3]), list(entries[3:])]
        c, d = gamma(e), gamma(cube(e))
        solutions = []
        for x in product(range(-2, 3), repeat=3):
            a = [list(x)] + e
            ac = cube(a)
            da, dc = det(a), det(ac)
            assert da == sum(ci*xi for ci, xi in zip(c, x))
            assert dc == sum(di*xi**3 for di, xi in zip(d, x))
            scaled = [[2*xi for xi in x]] + e
            for i in [1, 2]:
                for j in range(3):
                    assert cofactor(scaled, i, j) == 2*cofactor(a, i, j)
                    assert cofactor(cube(scaled), i, j) == 8*cofactor(ac, i, j)
            if da == dc == 1:
                simultaneous += 1
                solutions.append(x)
                for i in range(3):
                    p = [cofactor(a, i, j) for j in range(3)]
                    q = [cofactor(ac, i, j) for j in range(3)]
                    grad = [3*qj*a[i][j]**2 for j, qj in enumerate(q)]
                    dependent = all(p[j]*grad[k] == p[k]*grad[j]
                                    for j, k in combinations(range(3), 2))
                    identity = all(p[j] == a[i][j]**2*q[j] for j in range(3))
                    assert dependent == identity
            count += 1
        # A genuine replacement changes only the top row, preserving BOTH
        # determinant-one equations. No positivity assertion is made here.
        if len(solutions) >= 2:
            a = [list(solutions[0])] + e
            ap = [list(solutions[-1])] + e
            assert a != ap and a[1:] == ap[1:]
            assert det(a) == det(cube(a)) == det(ap) == det(cube(ap)) == 1
            replacements += 1
    assert count == 91125 and simultaneous > 0 and replacements > 0
    print('REPLACEMENT: all', count, 'small matrices satisfy Laplace and cofactor scaling')
    print('SINGULARITY: criterion checked on', simultaneous, 'simultaneous det-one matrices')
    print('REPLACEMENT: determinant-preserving changed-row pairs for', replacements, 'lower blocks')


def explicit_steps_and_border_tests():
    e = [[1,1,0], [1,1,1]]
    a = [[1,0,0]] + e
    ap = [[1,0,5]] + e
    assert det(a) == det(cube(a)) == det(ap) == det(cube(ap)) == 1
    print('ONE REPLACEMENT:', a, '->', ap, '; determinant pairs (1,1) -> (1,1)')
    ep, k1, k3 = border(e, e[0], [1,0], [1,0], [0,2], 1)
    assert k1 == k3 == 1
    extended = [ap[0] + [7]] + ep
    assert det(extended) == det(cube(extended)) == 1
    assert gamma(ep) == gamma(e) + [0]
    assert gamma(cube(ep)) == gamma(cube(e)) + [0]
    print('ONE UNIT BORDER:', extended, '; determinant pair (1,1), k1=k3=1')
    count = 0
    # Source Lemma 7.1, including arbitrary complement values and signs:
    # every E in {-1,0,1}^{2x3}, r=t E_i, t in {-2,-1,0,1,2},
    # w in {-1,0,1}^2, z in {-1,0,1}, both possible row indices i.
    for entries in product(range(-1,2), repeat=6):
        e = [list(entries[:3]), list(entries[3:])]
        for i, t, w, z in product(range(2), range(-2,3), product(range(-1,2), repeat=2), range(-1,2)):
            u, v = [0,0], [0,0]
            u[i], v[i] = t, t**3
            border(e, [t*a for a in e[i]], u, v, list(w), z)
            count += 1
    assert count == 196830
    print('FIDELITY: source Lemma 7.1 holds in all', count, 'specified exact bordering cases')
    print('SCOPE: finite algebraic checks; no positive order-ten witness or unboundedness proof')


if __name__ == '__main__':
    source_fidelity()
    all_small_replacement_identities()
    explicit_steps_and_border_tests()
    print('ALL_CHECKS_OK')
