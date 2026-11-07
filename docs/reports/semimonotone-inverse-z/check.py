#!/usr/bin/env python3
"""Exact certificates for Chauhan--Dubey Remark 4.1, second bullet.

Requires Python 3 and SymPy. No numerical eigensolver or floating-point test
is used. The finite checks complement the parameter identities and the
maximal-coordinate argument; they are not a kernel proof.
"""
from itertools import combinations, product
import sys
import sympy as sp


def principal(A, indices):
    return A.extract(indices, indices)


def source_fidelity():
    # Literal matrices from arXiv:2603.00639v1, Example 3.1.
    A = sp.Matrix([[0, -1, -1], [-2, 0, -1], [-3, -4, 0]])
    B = sp.Matrix([[1, -2, -2], [-2, 1, -2], [-2, -2, 1]])
    for label, M, strict in [('A', A, False), ('B', B, True)]:
        for i in range(3):
            assert M[i, i] > 0 if strict else M[i, i] >= 0
        count = 0
        for r in (2, 3):
            for indices in combinations(range(3), r):
                N = principal(M, indices)
                assert all(t < 0 for t in N * sp.ones(r, 1))
                count += 1
        print(f'FIDELITY Example 3.1 {label}: singleton certificates; '
              f'{count} larger principal blocks fail via x=1; '
              f'exact order 2 ({"strict" if strict else "ordinary"}); det={M.det()}')
    assert A.det() == -11 and B.det() == -27


def nested_product(x, y):
    """Coefficient product in the basis I,P,Q, with P^2=P,PQ=QP=Q,Q^2=Q."""
    a, b, c = x
    u, v, w = y
    return (a*u, a*v+b*u+b*v, a*w+c*u+b*w+c*v+c*w)


def parameter_identities():
    s, m = sp.symbols('s m', integer=True, positive=True)
    L = 3*s**2
    n = m*L
    d, a, b = 3*s*(s-1), 3*s-1, 3*s-2
    C = 1+n*b
    assert sp.expand(d-(s-1)*a) == s-1
    assert sp.expand(d-s*b) == -s
    # Actual diagonal, within-cluster and between-cluster entries.
    coeff = (L-1, -L, 1-C)
    inv_coeff = (1/(L-1), -1-1/(L-1), 1-1/C)
    assert sp.simplify(coeff[0]+coeff[1]/L+coeff[2]/n-d) == 0
    assert sp.simplify(coeff[1]/L+coeff[2]/n+a) == 0
    assert sp.simplify(coeff[2]/n+b) == 0
    assert tuple(map(sp.simplify, nested_product(coeff, inv_coeff))) == (1, 0, 0)
    assert tuple(map(sp.simplify, nested_product(inv_coeff, coeff))) == (1, 0, 0)
    assert sp.simplify(inv_coeff[2]/n-b/C) == 0
    print('SYMBOLIC: small-block margin=s-1; large-block row-sum bound=-s')
    print('SYMBOLIC: actual entries and both inverse products in the nested projection algebra verified')
    print('SYMBOLIC: between-cluster inverse entry=b/C')


def family_matrix(s, m):
    L = 3*s*s
    n = m*L
    d, a, b = 3*s*(s-1), 3*s-1, 3*s-2
    A = sp.Matrix(n, n, lambda i, j: d if i == j else (-a if i//L == j//L else -b))
    P = sp.Matrix(n, n, lambda i, j: sp.Rational(1, L) if i//L == j//L else 0)
    Q = sp.ones(n)/n
    return A, P, Q


def full_matrix_checks(s, m):
    L, n = 3*s*s, 3*m*s*s
    b = 3*s-2
    C = 1+n*b
    A, P, Q = family_matrix(s, m)
    I = sp.eye(n)
    assert P*P == P and Q*Q == Q and P*Q == Q and Q*P == Q
    E, F = I-P, P-Q
    assert E*E == E and F*F == F
    assert E*F == sp.zeros(n) and E*Q == sp.zeros(n) and F*Q == sp.zeros(n)
    assert (E.trace(), F.trace(), Q.trace()) == (n-m, m-1, 1)
    assert A == (L-1)*E-F-C*Q
    inverse_formula = E/(L-1)-F-Q/C
    # A separate exact inverse algorithm checks the proposed inverse.
    assert A.inv(method='DM') == inverse_formula
    assert A*inverse_formula == I and inverse_formula*A == I
    determinant = A.det(method='domain-ge')
    expected = (-1)**m*(L-1)**(n-m)*C
    assert determinant == expected
    assert inverse_formula[0, L] == sp.Rational(b, C) > 0
    assert determinant > 0 if m % 2 == 0 else determinant < 0
    if (s, m) == (2, 2):
        t = sp.Symbol('t')
        assert A.charpoly(t).as_expr().expand() == ((t-11)**22*(t+1)*(t+97)).expand()
        assert determinant == 97*11**22
        assert inverse_formula[0, L] == sp.Rational(4, 97)
    print(f'FULL s={s} m={m} n={n} k={n-s}: '
          f'det={(-1)**m}*{L-1}^{n-m}*{C}; '
          f'inertia=({n-m} positive,{m} negative,0 zero); '
          f'inverse_between={b}/{C}')


def occupancy_checks(s, m):
    """Exhaust all principal-block types, including ALL orders above the cutoff.

    A type records how many indices are chosen from each cluster. Permutation
    symmetry covers every actual principal subset, without sampling subsets.
    The small-block certificate uses a positive maximal coordinate of x>=0:
    (Bx)_i >= [d-(r-1)a] x_i > 0. For larger blocks B*1<0.
    """
    L = 3*s*s
    d, a, b = 3*s*(s-1), 3*s-1, 3*s-2
    small = large = 0
    for counts in product(range(L+1), repeat=m):
        r = sum(counts)
        if r == 0:
            continue
        if r <= s:
            assert d-(r-1)*a >= s-1 > 0
            small += 1
        else:
            bound = d-(r-1)*b
            assert bound <= -s < 0
            row_sums = [d-(c-1)*a-(r-c)*b for c in counts if c]
            assert all(t <= bound < 0 for t in row_sums)
            large += 1
    assert small+large == (L+1)**m-1
    print(f'PRINCIPAL TYPES s={s} m={m}: {small} nonempty small types; '
          f'{large} larger types; strict/ordinary cutoff={s}')


def main():
    print(f'Exact arithmetic: Python {sys.version_info.major}.{sys.version_info.minor}; SymPy {sp.__version__}')
    source_fidelity()
    parameter_identities()
    for s, m in ((2, 2), (2, 3), (3, 2)):
        full_matrix_checks(s, m)
        occupancy_checks(s, m)
    # Additional integer parameter checks, independent of matrix construction.
    for s in range(2, 21):
        for m in range(2, 21):
            L, n = 3*s*s, 3*m*s*s
            d, a, b = 3*s*(s-1), 3*s-1, 3*s-2
            C = 1+n*b
            assert n-s >= 1 and n-s < n and d > 0 and a > b > 0
            assert d-(s-1)*a == s-1 > 0 and d-s*b == -s < 0
            assert C > 1 and sp.Rational(b, C) > 0
    print('PARAMETERS: all 361 pairs 2<=s,m<=20 satisfy the exact sign and cutoff identities')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
