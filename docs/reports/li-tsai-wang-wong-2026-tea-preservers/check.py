#!/usr/bin/env python3
"""Exact finite checks for the TEA-preserver preregistration.

Gaussian integers are pairs of integers. Squared sup norms are integers.
For nonnegative squared norms A, B, C, sqrt(C) = sqrt(A) + sqrt(B)
iff C - A - B >= 0 and (C - A - B)**2 = 4*A*B.
No finite screen proves the arbitrary C*-algebra assertion.
"""
from fractions import Fraction
from itertools import product

ZERO = (0, 0)
ONE = (1, 0)
I = (0, 1)
UNITS_ZERO = (ZERO, ONE, (-1, 0), I, (0, -1))
GAUSSIAN = tuple(product(range(-1, 2), repeat=2))


def add(z, w):
    return (z[0] + w[0], z[1] + w[1])


def mul(z, w):
    return (z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0])


def conj(z):
    return (z[0], -z[1])


def abs2(z):
    return z[0]*z[0] + z[1]*z[1]


def vsum(x, y):
    return tuple(add(a, b) for a, b in zip(x, y))


def sup2(x):
    return max(map(abs2, x))


def tea2(a, b, c):
    delta = c - a - b
    return delta >= 0 and delta*delta == 4*a*b


def parallel(x, y):
    return any(abs2(a) == sup2(x) and abs2(b) == sup2(y)
               for a, b in zip(x, y))


def row_apply(row, x):
    return add(mul(row[0], x[0]), mul(row[1], x[1]))


def apply(rows, x):
    return tuple(row_apply(row, x) for row in rows)


def fidelity():
    # Source Lemma 2.1(a.2): beta=1/2. Use 2T to avoid denominators.
    rows = (((2, 0), ONE), (ONE, (2, 0)))
    vectors = tuple(product(GAUSSIAN, repeat=2))
    for x in vectors:
        tx = apply(rows, x)
        assert abs2(tx[0])-abs2(tx[1]) == 3*(abs2(x[0])-abs2(x[1]))
    for x, y in product(vectors, repeat=2):
        assert parallel(x, y) == parallel(apply(rows, x), apply(rows, y))
    # Complex phase distinguishes parallel from TEA, as the source emphasizes.
    x, y = (ONE, I), (ONE, (0, -1))
    assert tea2(sup2(x), sup2(y), sup2(vsum(x, y)))
    tx, ty = apply(rows, x), apply(rows, y)
    assert not tea2(sup2(tx), sup2(ty), sup2(vsum(tx, ty)))
    assert (sup2(tx), sup2(ty), sup2(vsum(tx, ty))) == (5, 5, 16)
    print('FIDELITY: Lemma 2.1(a.2), beta=1/2: 6561 parallel equivalences; phase TEA failure (5,5,16)')


def exact_screen(alphabet, row_count):
    vectors = tuple(product(UNITS_ZERO, repeat=2))
    pairs = [(x, y, vsum(x, y)) for x, y in product(vectors, repeat=2)
             if tea2(sup2(x), sup2(y), sup2(vsum(x, y)))]
    points = tuple(dict.fromkeys([x for triple in pairs for x in triple]
                                + list(product(GAUSSIAN, repeat=2))))
    point_index = {x: k for k, x in enumerate(points)}
    indices = [tuple(point_index[x] for x in triple) for triple in pairs]
    rows_all = tuple(product(alphabet, repeat=2))
    row_values = {row: tuple(abs2(row_apply(row, x)) for x in points)
                  for row in rows_all}
    count, survivors = 0, 0
    for rows in product(rows_all, repeat=row_count):
        count += 1
        values = tuple(max(v) for v in zip(*(row_values[row] for row in rows)))
        if not all(tea2(values[i], values[j], values[k]) for i, j, k in indices):
            continue
        survivors += 1
        # A genuinely independent structural check for this coefficient family:
        # all nonzero rows use one coordinate, both coordinates attain one scale.
        assert all(row[0] == ZERO or row[1] == ZERO for row in rows), rows
        scales = [max(abs2(row[j]) for row in rows) for j in range(2)]
        assert scales[0] == scales[1], rows
        scale2 = scales[0]
        assert all(values[k] == scale2*sup2(x) for k, x in enumerate(points)), rows
    print(f'SETTLEMENT_SCREEN: rows={row_count}, alphabet={len(alphabet)}, maps={count}, input_TEA_pairs={len(pairs)}, scalar_sup_seminorms={survivors}')
    return survivors


def phase_perturbation():
    # Lemma 1's hypothetical affine face slope, tested at rational parameters.
    alpha, beta = Fraction(1), Fraction(1, 4)
    s, t = Fraction(-16), Fraction(1, 16)
    assert s < -alpha/(2*beta*beta)
    x, y = (Fraction(1), t), (Fraction(0), s*beta*t)
    assert abs2(y) <= abs2(x)
    # p(x,y)=(alpha*|x|^2 + Re(beta*y*conj(x)))/|x|.
    numerator = alpha*abs2(x) + mul((beta, 0), mul(y, conj(x)))[0]
    p2 = numerator*numerator/abs2(x)
    assert numerator > 0
    assert p2 == Fraction(256, 257) < alpha*alpha
    # t -> -t has the same value, but their midpoint is ((1,0),(0,0)).
    xn, yn = (Fraction(1), -t), (Fraction(0), -s*beta*t)
    assert vsum((x, y), (xn, yn)) == ((2, 0), (0, 0))
    assert alpha*abs2(xn)+mul((beta, 0), mul(yn, conj(xn)))[0] == numerator
    print('LEMMA_1: rational phase perturbation p(v_t)^2=p(v_-t)^2=256/257 < p(midpoint)^2=1')


def nonunital_cutoff_model():
    # Exact finite-support c_0 model of the nine-corner boundedness argument.
    g = [tuple(Fraction(i == j) for i in range(4)) for j in range(3)]
    f = g
    e = [tuple(f[j][i]+(1-sum(f[k][i] for k in range(3)))/3
               for i in range(4)) for j in range(3)]
    assert all(sum(e[j][i] for j in range(3)) == 1 for i in range(4))
    a = (Fraction(2), Fraction(-3), Fraction(1), Fraction(5))
    corners = []
    for i, j in product(range(3), repeat=2):
        k = next(k for k in range(3) if k != i and k != j)
        x = tuple(e[i][r]*a[r]*e[j][r] for r in range(4))
        assert all(x[r]*g[k][r] == 0 for r in range(4))
        corners.append(x)
    assert tuple(sum(x[r] for x in corners) for r in range(4)) == a
    print('NONUNITAL_MODEL: rational local cutoffs sum to one; all nine orthogonal corners reconstruct a')


def main():
    fidelity()
    assert exact_screen(GAUSSIAN, 2) == 65
    assert exact_screen(UNITS_ZERO, 1) == 1
    exact_screen(UNITS_ZERO, 3)
    phase_perturbation()
    nonunital_cutoff_model()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
