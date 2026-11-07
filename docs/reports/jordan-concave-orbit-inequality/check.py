#!/usr/bin/env python3
"""Exact finite checks for Jeong's published Conjecture 4.2.

Python standard library only. No floating point, numerical tolerance, or random
sampling is used. These checks do not prove the unrestricted conjecture.
"""
from fractions import Fraction as Q
from itertools import product
from math import isqrt


def rational_sqrt(x):
    assert x >= 0
    a, b = isqrt(x.numerator), isqrt(x.denominator)
    return Q(a, b) if a*a == x.numerator and b*b == x.denominator else None


# Symmetric 2 by 2 matrix represented by (a11, a12, a22).
def eigenvalues(a):
    x, z, y = map(Q, a)
    gap = rational_sqrt((x-y)**2 + 4*z*z)
    return None if gap is None else ((x+y+gap)/2, (x+y-gap)/2)


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def abs_eigs(es):
    return tuple(sorted(map(abs, es), reverse=True))


FUNCTIONS = [
    ('identity', lambda t: t),
    ('constant_one', lambda t: Q(1)),
    ('resolvent', lambda t: t/(1+t)),
    ('one_plus_resolvent', lambda t: 1+t/(1+t)),
    ('min_one', lambda t: min(t, Q(1))),
    ('min_two', lambda t: min(t, Q(2))),
    ('positive_indicator', lambda t: Q(t > 0)),
    ('one_plus_positive_indicator', lambda t: 1+Q(t > 0)),
]


def certify_two_orbits(alpha, beta, target):
    """An exact algebraic witness for 2 by 2 orthogonal-orbit domination.

    C = diag(alpha1,alpha2), D = R diag(beta1,beta2) R^T,
    R has cos^2(theta)=t. Compute det(C+D) with the square of the
    off-diagonal entry, avoiding any irrational square roots. An orthogonal
    diagonalization then aligns C+D with the target. Checking both ordered
    eigenvalue gaps certifies cone domination after this common alignment.
    """
    a1, a2 = alpha
    b1, b2 = beta
    c1, c2 = target
    assert a1 >= a2 >= 0 and b1 >= b2 >= 0 and c1 >= c2 >= 0
    total = a1+a2+b1+b2
    lower = max(a1+b2, a2+b1)
    upper = a1+b1
    z = max(lower, c1)
    assert z <= min(upper, total-c2), (alpha, beta, target)
    determinant = z*(total-z)
    denominator = (a1-a2)*(b1-b2)
    det_at_zero = (a1+b2)*(a2+b1)
    t = (det_at_zero-determinant)/denominator if denominator else Q(0)
    assert 0 <= t <= 1
    d11, d22 = b1*t+b2*(1-t), b1*(1-t)+b2*t
    off_diagonal_squared = (b1-b2)**2*t*(1-t)
    assert (a1+d11)*(a2+d22)-off_diagonal_squared == determinant
    assert a1+a2+d11+d22 == total
    assert z >= total-z and z-c1 >= 0 and total-z-c2 >= 0


def fidelity():
    # Jeong Theorem 3.5, equality x=e1, y=e2, p=1/2.
    # ||x||_{1/2} = (sum sqrt(abs(eigenvalue)))^2.
    def schatten_half(es):
        roots = [rational_sqrt(abs(Q(t))) for t in es]
        assert all(t is not None for t in roots)
        return sum(roots, Q(0))**2
    x, y, xy = (1, 0), (0, 1), (1, 1)
    lhs = schatten_half(xy)
    rhs = 2*(schatten_half(x)+schatten_half(y))
    assert lhs == rhs == 4
    print('FIDELITY Theorem 3.5: p=1/2, x=e1, y=e2, lhs=rhs=4')


def diagonal_checks():
    count = 0
    for rank in (1, 2, 3):
        vectors = list(product(range(-2, 3), repeat=rank))
        for a, b in product(vectors, repeat=2):
            for name, f in FUNCTIONS:
                # The identity automorphisms witness the full cone inequality
                # in the componentwise Jordan algebra, with no permutations.
                assert all(f(Q(abs(x+y))) <= f(Q(abs(x)))+f(Q(abs(y)))
                           for x, y in zip(a, b)), (rank, a, b, name)
                count += 1
    print('DIAGONAL ranks=1..3, entries=-2..2, functions=8, checks='+str(count))


def matrix_checks():
    matrices = [tuple(map(Q, a)) for a in product(range(-2, 3), repeat=3)]
    matrices = [(a, eigenvalues(a)) for a in matrices if eigenvalues(a) is not None]
    pairs = noncommuting = singular_sums = checks = excluded = 0
    for (a, ea), (b, eb) in product(matrices, repeat=2):
        ec = eigenvalues(add(a, b))
        if ec is None:
            excluded += 1
            continue
        pairs += 1
        # [A,B]12 = (a11-a22)b12 - (b11-b22)a12.
        noncommuting += (a[0]-a[2])*b[1] != (b[0]-b[2])*a[1]
        singular_sums += 0 in ec
        for name, f in FUNCTIONS:
            alpha = tuple(map(f, abs_eigs(ea)))
            beta = tuple(map(f, abs_eigs(eb)))
            target = tuple(map(f, abs_eigs(ec)))
            certify_two_orbits(alpha, beta, target)
            # Jeong Theorem 3.4, checked independently as a trace scalar.
            assert sum(target) <= sum(alpha)+sum(beta)
            checks += 1
    assert noncommuting > 0 and singular_sums > 0
    print('MATRIX entries=-2..2, rational-spectrum matrices='+str(len(matrices)))
    print('MATRIX pairs='+str(pairs)+', noncommuting='+str(noncommuting)+
          ', singular_sums='+str(singular_sums)+', checks='+str(checks))
    print('EXCLUDED pairs with irrational eigenvalues of A+B='+str(excluded))


def contraction_boundary():
    # Scalar instances of the Jensen step, including contractions of norm <1,
    # a discontinuity at zero and f(0)>0; all values are rational.
    count = 0
    grid = [Q(i, 4) for i in range(9)]
    for s in grid:
        for t in [Q(i, 8) for i in range(9)]:
            for name, f in FUNCTIONS:
                assert t*f(s) <= f(t*s)-(1-t)*f(Q(0)), (s, t, name)
                count += 1
    print('JENSEN scalar contraction checks='+str(count))
    print('BOUNDARIES include zero, singular sums, f(0)>0, endpoint-discontinuous f')


def main():
    fidelity()
    diagonal_checks()
    matrix_checks()
    contraction_boundary()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
