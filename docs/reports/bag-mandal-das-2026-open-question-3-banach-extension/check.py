#!/usr/bin/env python3
"""Exact finite checks for the Bag--Mandal--Das Banach-extension draft.

These checks test source fidelity and identities used in the proposed general
argument. They do not prove a statement about all Banach spaces. All arithmetic
used as evidence is integral or rational; the norm checks use squared Euclidean
norms or the rational operator norm induced by the maximum norm.
"""
from fractions import Fraction as Q
from itertools import product


def eye(n):
    return tuple(tuple(Q(i == j) for j in range(n)) for i in range(n))


def zero(n):
    return tuple(tuple(Q(0) for _ in range(n)) for _ in range(n))


def add(a, b):
    return tuple(tuple(x + y for x, y in zip(ar, br))
                 for ar, br in zip(a, b))


def scale(c, a):
    return tuple(tuple(c * x for x in row) for row in a)


def sub(a, b):
    return add(a, scale(-1, b))


def mul(a, b):
    n = len(a)
    return tuple(tuple(sum((a[i][k] * b[k][j] for k in range(n)), Q(0))
                       for j in range(n)) for i in range(n))


def power(a, n):
    assert n >= 0
    out = eye(len(a))
    for _ in range(n):
        out = mul(out, a)
    return out


def apply(a, x):
    return tuple(sum((v * w for v, w in zip(row, x)), Q(0)) for row in a)


def determinant2(a):
    assert len(a) == 2
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def inverse2(a):
    det = determinant2(a)
    assert det != 0
    return scale(1 / det, ((a[1][1], -a[0][1]), (-a[1][0], a[0][0])))


def max_norm(a):
    return max((sum((abs(x) for x in row), Q(0)) for row in a), default=Q(0))


def all_two_matrices():
    for v in product((Q(-1), Q(0), Q(1)), repeat=4):
        yield ((v[0], v[1]), (v[2], v[3]))


def check_source_fidelity():
    # Source example labelled finitedim_example1: T = diag(0, 2, 3).
    a = ((Q(0), Q(0), Q(0)), (Q(0), Q(2), Q(0)),
         (Q(0), Q(0), Q(3)))
    count = 0
    for n in range(1, 13):
        an = power(a, n)
        for x in product(range(-2, 3), repeat=3):
            y = apply(an, x)
            assert y == (0, 2 ** n * x[1], 3 ** n * x[2])
            assert sum(t * t for t in y) == 4 ** n * x[1] ** 2 + 9 ** n * x[2] ** 2
            count += 1
    assert apply(a, (1, 0, 0)) == (0, 0, 0)
    print(f'SOURCE_FIDELITY diag(0,2,3): {count} exact power/norm-square cases OK')


def check_difference_of_powers():
    scalars = (Q(-2), Q(-1), Q(-1, 2), Q(0), Q(1, 2), Q(1), Q(2))
    mats = [()] + [((s,),) for s in scalars] + list(all_two_matrices())
    mats += [eye(3), zero(3), ((Q(0), Q(1), Q(0)), (Q(0), Q(0), Q(1)),
                             (Q(0), Q(0), Q(0))),
             ((Q(1), Q(2), Q(3)), (Q(-1), Q(0), Q(2)),
              (Q(2), Q(-2), Q(1)))]
    count = 0
    for a in mats:
        for n in range(1, 7):
            for lam in scalars:
                s = zero(len(a))
                for k in range(n):
                    s = add(s, scale(lam ** (n - 1 - k), power(a, k)))
                factor = sub(scale(lam, eye(len(a))), a)
                lhs = sub(scale(lam ** n, eye(len(a))), power(a, n))
                assert mul(factor, s) == lhs == mul(s, factor)
                count += 1
    print(f'DIFFERENCE_OF_POWERS: {count} exact cases (dimensions 0..3, N=1..6) OK')


def check_resolvent_branches():
    # Small matrices B with ||B||_infinity <= 1/4. Use c=1/2.
    # For nonsingular B, T=B^-1 has ||T^-1|| <= 1/4; use d=2.
    outside = (Q(-2), Q(-1), Q(-1, 2), Q(1, 2), Q(1), Q(2))
    inside = (Q(-1), Q(-1, 2), Q(0), Q(1, 2), Q(1))
    contraction_cases = expansion_cases = 0
    for m in all_two_matrices():
        b = scale(Q(1, 8), m)
        assert max_norm(b) <= Q(1, 4)
        for n in range(1, 7):
            bn = power(b, n)
            assert max_norm(bn) <= Q(1, 4) ** n
            for lam in outside:
                assert abs(lam) > Q(1, 4)
                assert determinant2(sub(scale(lam ** n, eye(2)), bn)) != 0
                contraction_cases += 1
        if determinant2(b) == 0:
            continue
        t = inverse2(b)
        assert mul(t, b) == eye(2) == mul(b, t)
        for n in range(1, 7):
            tn = power(t, n)
            assert mul(tn, power(b, n)) == eye(2)
            for lam in inside:
                assert abs(lam) < 2
                assert determinant2(sub(scale(lam ** n, eye(2)), tn)) != 0
                expansion_cases += 1
    print(f'RESOLVENT_BRANCHES: {contraction_cases} contracting and {expansion_cases} expanding exact rational cases OK')


def check_mixed_gap_obstruction():
    # For diag(1/2,2), squared coordinate weights p,1-p can be chosen
    # exactly so that both the input and its n-th image have Euclidean norm 1.
    for n in range(1, 33):
        small, large = Q(1, 4) ** n, Q(4) ** n
        second_weight = (1 - small) / (large - small)
        assert 0 < second_weight < 1
        assert small * (1 - second_weight) + large * second_weight == 1
    # This diagnoses the separate converse assertion, not Question 3.
    print('MIXED_GAP_OBSTRUCTION diag(1/2,2): 32 exact squared-weight cases OK')


def main():
    check_source_fidelity()
    check_difference_of_powers()
    check_resolvent_branches()
    check_mixed_gap_obstruction()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
