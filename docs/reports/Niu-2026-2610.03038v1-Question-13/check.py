#!/usr/bin/env python3
"""Exact, independent fidelity and preregistration checks; no Lean proof.

Source: Tong Niu, arXiv:2610.03038v1, Equations (1)--(3), Question 13.
Run: python3 check.py. Requires SymPy. All evidence uses exact arithmetic.
"""
from fractions import Fraction
from math import factorial, comb, prod
import json
import sympy as sp


def recurrence_values(limit):
    values = {1: 0, 2: 2}
    for n in range(3, limit + 1):
        numerator = ((n*n - 1)*(364*n*n - 503*n + 153)*values[n-1]
                     + 3*(3*n-1)*(3*n-2)*(2*n-3)*(7*n-2)*values[n-2])
        denominator = 2*(n+1)**2*(n+2)*(7*n-9)
        assert denominator > 0 and numerator % denominator == 0
        values[n] = numerator // denominator
    return values


def enumerate_tableaux(n):
    """Literal enumeration by adjoining the next label to an allowed cell."""
    rows = [[] for _ in range(3)]
    total = odd = 0
    def visit(label):
        nonlocal total, odd
        if label == 3*n + 1:
            assert sorted(x for row in rows for x in row) == list(range(1, 3*n+1))
            assert all(rows[i][j] < rows[i+1][j]
                       for i in range(2) for j in range(n))
            total += 1
            odd += rows[1][0] % 2
            return
        for row in range(3):
            if len(rows[row]) == n or (row and len(rows[row]) >= len(rows[row-1])):
                continue
            rows[row].append(label)
            visit(label+1)
            rows[row].pop()
    visit(1)
    return odd, total


def finite_sum(n):
    return sum((Fraction(8*m*(m+1)*factorial(3*n-2*m-1),
                         factorial(n)*factorial(n+2)*factorial(n-2*m))
                for m in range(1, n//2+1)), Fraction(0))


def main():
    print('Python exact arithmetic; SymPy', sp.__version__)
    a = recurrence_values(140)
    print('a_n from Theorem 2 for 1 <= n <= 40:')
    print(json.dumps({n:a[n] for n in range(1,41)},sort_keys=True))
    assert [a[n] for n in range(1,11)] == [0,2,16,182,2400,35310,562848,9540674,169777504,3142665968]
    for n in range(1,6):
        odd,total = enumerate_tableaux(n)
        assert odd == a[n]
        assert total == 2*factorial(3*n)//(factorial(n)*factorial(n+1)*factorial(n+2))
        print(f'BRUTE n={n}: odd={odd}, total={total}, recurrence_matches=True')
    assert all(finite_sum(n) == a[n] for n in range(1,141))
    print('SOURCE_FINITE_SUM_MATCH 1 <= n <= 140')

    # Exact linear algebra over Q: full column rank excludes complex coefficients too.
    # A finite window does not exclude a relation that begins beyond this window.
    for start, stop in [(2,39),(20,100)]:
        readings=[]
        for degree in range(13):
            rows = [[sp.Integer(n)**k * a[n+shift]
                     for shift in (1,0) for k in range(degree+1)]
                    for n in range(start,stop+1)]
            mat=sp.polys.matrices.DomainMatrix.from_Matrix(sp.Matrix(rows))
            rank=mat.rank()
            assert rank == 2*(degree+1)
            readings.append({'degree':degree,'rank':rank,'columns':2*(degree+1),'nullity':0})
        print(f'FIRST_ORDER_FIT n={start}..{stop}:',json.dumps(readings))
    print('FINITE_WINDOW_LIMIT: coefficient degrees <= 12 only; no unbounded or eventual conclusion from fits')

    n,m=sp.symbols('n m')
    c0=-3*(3*n+5)*(3*n+4)*(2*n+1)*(7*n+12)
    c1=-(n+1)*(n+3)*(364*n*n+953*n+603)
    c2=2*(n+3)**2*(n+4)*(7*n+5)
    r=-(2*n+1)*(8*n+15)/(2*(n+3)*(8*n+7))
    residual=sp.cancel(c0+c1*r+c2*r*r.subs(n,n+1))
    assert residual == 0
    print('HYPERGEOMETRIC_SOLUTION_RESIDUAL',residual)
    def b(k):
        return Fraction((-1)**k*(8*k+7)*comb(2*k,k),4**k*(k+1)*(k+2))
    assert all(b(k) != 0 and b(k+1)/b(k) == r.subs(n,k) for k in range(141))
    assert all(a[k] > 0 and r.subs(n,k) < 0 for k in range(2,140))
    assert all(Fraction(a[k+1],1)*b(k) != Fraction(a[k],1)*b(k+1) for k in range(2,140))
    w=-6*(3*n+4)*(3*n+5)*(7*n+12)*(8*n+7)/((n+3)*(2*n+3)*(7*n+5)*(8*n+23))
    residual=sp.cancel(w-c0/(c2*r*r.subs(n,n+1)))
    assert residual == 0
    print('ORDER_REDUCTION_RESIDUAL',residual)
    A=-6*(3*n+4)*(3*n+5)*(8*n+7)
    B=(n+3)*(2*n+3)*(8*n+23)
    C=7*n+5
    residual=sp.cancel(w-(A/B)*C.subs(n,n+1)/C)
    assert residual == 0
    print('GOSPER_FACTORIZATION_RESIDUAL',residual)
    ra=[-sp.Rational(4,3),-sp.Rational(5,3),-sp.Rational(7,8)]
    rb=[-sp.Integer(3),-sp.Rational(3,2),-sp.Rational(23,8)]
    assert all(A.subs(n,t)==0 for t in ra) and all(B.subs(n,t)==0 for t in rb)
    gaps=[v-u for u in ra for v in rb]
    assert all(not (gap.is_integer and gap>=0) for gap in gaps)
    print('ROOT_B_MINUS_ROOT_A',[str(t) for t in gaps])
    lc=sp.Poly(A,n).LC()-sp.Poly(B,n).LC()
    assert lc == -448
    print('LEADING_COEFFICIENT_DIFFERENCE',lc,'; polynomial degree obstruction d+3 > 1')
    # This identity is the source certificate; independently transcribed, not executed from source.
    p0=-2*(n+1)**2*(n+2)*(7*n-9)
    p1=(n*n-1)*(364*n*n-503*n+153)
    p2=3*(3*n-1)*(3*n-2)*(2*n-3)*(7*n-2)
    t=3*n-2*m-1
    u1=n*(n+2)*(n-2*m)/prod(t-i for i in range(3))
    u2=n*(n-1)*(n+1)*(n+2)*(n-2*m)*(n-2*m-1)/prod(t-i for i in range(6))
    v=(m+2)*(n-2*m)*(n-2*m-1)/(m*t*(t-1))
    W=(4*(7*n-9)*m**3-2*(7*n-9)*(16*n-11)*m**2
       +(539*n**3-1526*n*n+1339*n-360)*m
       -3*(n-1)*(147*n**3-392*n*n+308*n-72))
    R=-2*(m-1)*(2*m+1)*(n+1)*(n+2)*W/prod(t-i for i in range(4))
    residual=sp.cancel(p0+p1*u1+p2*u2-R.subs(m,m+1)*v+R)
    assert residual==0
    print('SOURCE_TELESCOPING_CERTIFICATE_RESIDUAL',residual)
    p,q,z=sp.symbols('p q z')
    total=p+q+z
    corner_ratios=((p+2)*(p-q)*(p-z+1)/(total*(p-q+1)*(p-z+2))
                   +(q+1)*(p-q+2)*(q-z)/(total*(p-q+1)*(q-z+1))
                   +z*(p-z+3)*(q-z+2)/(total*(p-z+2)*(q-z+1)))
    residual=sp.cancel(corner_ratios-1)
    assert residual==0
    print('THREE_ROW_CORNER_INDUCTION_RESIDUAL',residual)
    print('SCOPE: exact finite counts and rational identities; neither claim nor its negation has been proved in Lean')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
