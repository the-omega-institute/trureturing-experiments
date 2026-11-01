#!/usr/bin/env python3
"""Exact fidelity and algebra checks for Doust--Weston Question 4.9.

Requires Python 3 and SymPy. No floating-point values or approximate signs are
used. These checks are ordinary algebraic evidence, not a Lean proof.
"""
from fractions import Fraction
from itertools import product
import sympy as S


def distance_matrix(A, B):
    return S.Matrix([[0, 1, A, B], [1, 0, 1, A],
                     [A, 1, 0, 1], [B, A, 1, 0]])


def quadratic(M, x):
    return (x.T * M * x)[0]


def fidelity():
    # Source Section 4.3: exact odd-cycle determinant and inverse-entry sum.
    for n in range(3, 24, 2):
        D = S.Matrix(n, n, lambda i, j: min(abs(i-j), n-abs(i-j)))
        assert D.det() == S.Rational(n*n-1, 4)
        one = S.ones(n, 1)
        assert (one.T * D.inv() * one)[0] == S.Rational(4*n, n*n-1)
    print('FIDELITY: 11 odd cycles C_3 through C_23; both source identities exact')

    # Source's C5 threshold: exact factorization, with the first vanishing
    # nonconstant eigenvalue at A = (3 + sqrt(5))/2, where A = 2**p.
    A, z = S.symbols('A z')
    D = S.Matrix(5, 5, lambda i, j: 0 if i == j else
                 (1 if min(abs(i-j), 5-abs(i-j)) == 1 else A))
    assert S.expand(D.det() - 2*(A+1)*(A*A-3*A+1)**2) == 0
    c = (3 + S.sqrt(5))/2
    assert S.simplify(c*c-3*c+1) == 0
    assert S.expand(D.charpoly(z).as_expr() -
                    (z-2*(A+1))*(z*z+(A+1)*z-(A*A-3*A+1))**2) == 0
    print('FIDELITY: source C5 threshold polynomial and full characteristic polynomial exact')


def candidate_symbolic():
    A, B, s, u, v = S.symbols('A B s u v', real=True)
    x = S.Matrix([s+u, -s+v, -s-v, s-u])
    Q = quadratic(distance_matrix(A, B), x)
    split = 2*(B-1-2*A)*s*s - 2*(B*u*u+2*(A-1)*u*v+v*v)
    assert sum(x) == 0
    assert S.expand(Q-split) == 0

    # Coordinates span the entire zero-sum subspace, not just test vectors.
    x0, x1, x2 = S.symbols('x0 x1 x2', real=True)
    arbitrary = S.Matrix([x0, x1, x2, -x0-x1-x2])
    coord = {s: -(x1+x2)/2, u: x0+(x1+x2)/2, v: (x1-x2)/2}
    assert (x.subs(coord)-arbitrary) == S.zeros(4, 1)
    endpoint = 2*A*(A-4)*s*s - 2*((A-1)*u+v)**2
    assert S.expand(Q.subs(B, (A-1)**2)-endpoint) == 0
    witness = S.Matrix([1, -(A-1), A-1, -1])
    assert sum(witness) == 0
    assert S.expand(quadratic(distance_matrix(A, B), witness)
                    - 2*((A-1)**2-B)) == 0
    # The cross block vanishes by reversal symmetry.
    e = S.Matrix([1, -1, -1, 1])
    odd = S.Matrix([u, v, -v, -u])
    assert S.expand((e.T*distance_matrix(A, B)*odd)[0]) == 0
    print('CANDIDATE: symbolic zero-sum decomposition, endpoint square and obstruction exact')

    # The derivative numerator after log(2**q-1) = q*log(2)
    # + log(1-2**(-q)); t denotes 2**(-q), 0<t<1.
    q, L, t, ell = S.symbols('q L t ell', positive=True)
    raw = 2*(q*L/(1-t) - (q*L+ell))/(q*q)
    positive_form = 2*(q*L*t/(1-t)-ell)/(q*q)
    assert S.factor(raw-positive_form) == 0
    print('CANDIDATE: derivative rewrite exact; analytic sign uses q>0 and 0<t<1')


def candidate_small_cases():
    # Exact geometric distances in units of a = 2*pi/(b+3).
    count = 0
    for denominator in range(1, 10):
        for numerator in range(denominator+1, 3*denominator):
            b = Fraction(numerator, denominator)
            circumference = b+3
            assert Fraction(1, 3) < Fraction(2, 1)/circumference < Fraction(1, 2)
            assert 0 < 3 < circumference
            expected = [[0, 1, 2, b], [1, 0, 1, 2],
                        [2, 1, 0, 1], [b, 2, 1, 0]]
            for i, j in product(range(4), repeat=2):
                assert min(abs(i-j), circumference-abs(i-j)) == expected[i][j]
            count += 1
    print(f'CANDIDATE: {count} rational gap parameters; all 16 geodesic distances exact')

    # Every A on this rational grid lies in (2,4). At B=(A-1)^2
    # the square certificate must be nonpositive for every grid vector.
    count = 0
    for k in range(1, 16):
        A = S.Rational(2) + S.Rational(k, 8)
        B = (A-1)**2
        M = distance_matrix(A, B)
        for ss, uu, vv in product(range(-2, 3), repeat=3):
            x = S.Matrix([ss+uu, -ss+vv, -ss-vv, ss-uu])
            value = quadratic(M, x)
            assert value <= 0
            count += 1
        witness = S.Matrix([1, -(A-1), A-1, -1])
        for ratio in [S.Rational(1, 4), S.Rational(1, 2), S.Rational(3, 4)]:
            assert quadratic(distance_matrix(A, ratio*B), witness) > 0
    print(f'CANDIDATE: {count} endpoint vectors and 45 strict obstruction cases exact')
    # A wrong sign in the endpoint square must be detected.
    assert quadratic(distance_matrix(3, 4), S.Matrix([0, 1, -1, 0])) == -2
    assert quadratic(distance_matrix(3, 4), S.Matrix([0, 1, -1, 0])) != 2
    print('NEGATIVE CONTROL: reversing the square sign is rejected')


if __name__ == '__main__':
    fidelity()
    candidate_symbolic()
    candidate_small_cases()
    print('ALL_CHECKS_OK')
