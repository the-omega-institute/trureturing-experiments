#!/usr/bin/env python3
"""Exact checks for Fisk Question 10 (0808.1850v1), not a Lean proof.

Fidelity: the repeated-factor example in math/0612833v2, Section 21.10.
Candidate: a genuine 5 by 5 coefficient-Toeplitz minor for every integer t >= 4.
Only exact SymPy algebra and integer arithmetic are used.
"""
from itertools import combinations, permutations
from math import comb
import sympy as sp

x, y, t = sp.symbols("x y t", real=True)
m = sp.symbols("m", nonnegative=True, integer=True)


def zero(expr):
    assert sp.cancel(sp.expand(expr)) == 0, expr


def determinant_by_permutations(matrix):
    n = matrix.rows
    result = 0
    for perm in permutations(range(n)):
        inversions = sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        result += (-1) ** inversions * sp.prod(matrix[i, perm[i]] for i in range(n))
    return sp.expand(result)


def minor_from_coefficients(f, rows, columns):
    coefficients = sp.Poly(f, y)
    return sp.Matrix([
        [coefficients.nth(c - r) if c >= r else 0 for c in columns]
        for r in rows
    ])


def fidelity():
    # The book displays phi((1+x+y)^4), whose minors are constants times powers
    # of (1+x). Enumerate every minor of its first 4 rows and first 7 columns.
    f = sp.expand((1 + x + y) ** 4)
    displayed = minor_from_coefficients(f, range(4), range(7))
    assert list(displayed.row(0)) == [sp.expand((1+x)**4),
                                    sp.expand(4*(1+x)**3), sp.expand(6*(1+x)**2),
                                    4*(1+x), 1, 0, 0]
    count = nonzero = 0
    for n in range(5):
        for rows in combinations(range(4), n):
            for columns in combinations(range(7), n):
                mat = displayed.extract(rows, columns)
                actual = determinant_by_permutations(mat)
                constants = sp.Matrix([
                    [comb(4, c-r) if 0 <= c-r <= 4 else 0 for c in columns]
                    for r in rows
                ])
                constant = determinant_by_permutations(constants)
                assert constant >= 0
                exponent = 4*n + sum(rows) - sum(columns)
                if constant:
                    assert exponent >= 0
                    zero(actual - constant*(1+x)**exponent)
                    nonzero += 1
                else:
                    zero(actual)
                count += 1
    assert count == 330
    print(f"FIDELITY: all {count} minors of the displayed 4x7 block checked; {nonzero} nonzero")
    print("FIDELITY: each nonzero minor is c*(1+x)^k with integer c>0 and k>=0; its roots are -1")


def candidate():
    f = sp.expand((t*x + y + 1)*(x + y + t))
    expected = t*x**2 + (t+1)*x*y + (t*t+1)*x + y*y + (t+1)*y + t
    zero(f - expected)
    coefficients = sp.Poly(f, x, y).coeffs()
    assert len(coefficients) == 6
    # A polynomial in m=t-4 with nonnegative coefficients and positive constant
    # is positive for all m>=0. This certifies the six coefficient signs.
    for coefficient in coefficients:
        p = sp.Poly(coefficient.subs(t, m+4), m)
        assert all(c >= 0 for c in p.all_coeffs()) and p.nth(0) > 0
    xr, xi, yr, yi = sp.symbols("xr xi yr yi", real=True)
    for affine, expected_im in [(t*x+y+1, t*xi+yi), (x+y+t, xi+yi)]:
        zero(sp.im(sp.expand(affine.subs({x:xr+sp.I*xi, y:yr+sp.I*yi}))) - expected_im)
    # Thus both imaginary parts are positive when xi,yi>0 and t>=4;
    # neither affine factor vanishes. No sample of complex points is used.
    u, v = t*x + 1, x + t
    a, b = sp.expand(u*v), sp.expand(u+v)
    rows, columns = list(range(5)), list(range(1, 6))
    mat = minor_from_coefficients(f, rows, columns)
    expected_mat = sp.Matrix([[b,1,0,0,0],[a,b,1,0,0],[0,a,b,1,0],
                              [0,0,a,b,1],[0,0,0,a,b]])
    assert mat == expected_mat
    actual = determinant_by_permutations(mat)
    zero(actual - mat.det(method="domain-ge"))
    recurrence = [sp.Integer(1), b]
    for n in range(2, 6):
        recurrence.append(sp.expand(b*recurrence[-1] - a*recurrence[-2]))
    for n in range(1, 6):
        zero(determinant_by_permutations(mat[:n,:n]) - recurrence[n])
    factorized = b*(b*b-a)*(b*b-3*a)
    zero(actual - factorized)
    zero(actual - (u+v)*(u*u+u*v+v*v)*(u*u-u*v+v*v))
    A, C, S = t*t-t+1, t*t-4*t+1, t*t-1
    Q = A*x*x - C*x + A
    zero(Q - (u*u-u*v+v*v))
    zero(4*A*A-C*C-3*S*S)
    for quantity in [A,C,S]:
        p = sp.Poly(quantity.subs(t,m+4),m)
        assert all(c >= 0 for c in p.all_coeffs()) and p.nth(0) > 0
    z = (C + sp.I*sp.sqrt(3)*S)/(2*A)
    zero(Q.subs(x,z))
    zero(sp.re(z) - C/(2*A))
    zero(actual.subs(x,z))
    print("FAMILY: exact coefficient positivity and upper-half-plane affine-factor identities checked")
    print("FAMILY: rows [0,1,2,3,4], columns [1,2,3,4,5]; permutation determinant equals continuant and factors")
    print("FAMILY: Q=(t^2-t+1)*x^2-(t^2-4*t+1)*x+(t^2-t+1)")
    print("FAMILY: Q(z_t)=det(M_t)(z_t)=0 identically; Re(z_t)=C/(2*A)>0 for all real t>=4")
    for value in range(4, 41):
        root = z.subs(t,value)
        zero(Q.subs({t:value,x:root}))
        zero(actual.subs({t:value,x:root}))
        assert sp.re(root) > 0
    print("SMALL CASES: all integers 4<=t<=40 checked again with exact algebraic roots")
    small = 5*(x+1)*(21*x*x+33*x+21)*(13*x*x-x+13)
    zero(actual.subs(t,4)-small)
    zero(small - (1365*x**5+3405*x**4+4605*x**3+4605*x**2+3405*x+1365))
    root4 = (1+15*sp.I*sp.sqrt(3))/26
    zero(z.subs(t,4)-root4)
    zero((13*x*x-x+13).subs(x,root4))
    assert sp.re(root4) == sp.Rational(1,26)
    zero((13*x*x-x+13).subs(x,x+sp.Rational(1,26)) - 13*(x*x+sp.Rational(675,676)))
    print("t=4: det=5*(x+1)*(21*x^2+33*x+21)*(13*x^2-x+13)")
    print("t=4: z=(1+15*I*sqrt(3))/26, Re(z)=1/26, exact determinant evaluation zero")
    print("t=4: translation to w^2+675/676 checked")


if __name__ == "__main__":
    print("SymPy", sp.__version__, "; exact arithmetic only")
    fidelity()
    candidate()
    print("ALL_CHECKS_OK")
