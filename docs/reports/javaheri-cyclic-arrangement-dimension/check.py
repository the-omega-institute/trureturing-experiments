#!/usr/bin/env python3
"""Exact checks for arXiv:2610.33653v1, Conjecture 1.3.

All evidence uses integers or exact symbolic algebra. This is independent
computational evidence for a preregistration, not a Lean proof.
"""
from fractions import Fraction
from itertools import product
import sympy as sp


def sides(xs, a, b):
    """Literal cyclic inequality, specialized only when c is a natural number."""
    c = Fraction(a + b + 1, 2)
    assert c.denominator == 1 and len(xs) >= 1
    assert all(x >= 0 for x in xs)
    cyc = sum(x ** a * xs[(i + 1) % len(xs)] ** b for i, x in enumerate(xs))
    lhs = sum(xs) * cyc
    rhs = sum(x ** c.numerator for x in xs) ** 2
    return lhs, rhs, cyc


def identity(lhs, rhs):
    assert sp.expand(lhs - rhs) == 0


def main():
    x, y, n, m, r, t = sp.symbols('x y n m r t')

    # Fidelity: Corollary 4.2(ii), label cor:examples, says N(0,b)=1 for b != 1.
    # Reproduce b=3. Every n>=2 has the witness (2,1,...,1), with D_n=2(n-1).
    assert sides([2], 0, 3)[:2] == (16, 16)
    identity((n + 1) * (n + 7) - (n + 3) ** 2, 2 * (n - 1))
    for dim in range(2, 101):
        lhs, rhs, _ = sides([2] + [1] * (dim - 1), 0, 3)
        assert lhs - rhs == 2 * (dim - 1) > 0
    print('FIDELITY cor:examples(ii): N(0,3)=1; D_n(2,1,...,1)=2(n-1), checked n=2..100')
    lhs, rhs, cyc = sides([1, 2], 0, 3)
    assert (lhs, rhs, cyc) == (27, 25, 9)
    print(f'FIDELITY n=2, x=(1,2), a=0,b=3,c=2: LHS={lhs}, RHS={rhs}, D={lhs-rhs}')

    # Global validity in dimension 2: verify the complete algebraic identity.
    gap = (x ** 4 + y ** 4) ** 2 - (x + y) * (x ** 6 * y + y ** 6 * x)
    factor = (x-y) ** 2 * (x*x+x*y+y*y) * (x ** 4-x*x*y*y+y ** 4)
    identity(gap, factor)
    identity(x ** 4-x*x*y*y+y ** 4, (x*x-y*y) ** 2+x*x*y*y)
    print('N(6,1)>=2: exact two-variable factorization and square decomposition verified')
    for u, v in product(range(21), repeat=2):
        lhs, rhs, _ = sides([u, v], 6, 1)
        assert lhs <= rhs
    for u in range(21):
        lhs, rhs, _ = sides([u], 6, 1)
        assert lhs == rhs
    print('ADDITIONAL finite checks: 441 nonnegative integer pairs and 21 singletons passed')

    # All n>=3: compute source-oriented cyclic sums, not a symmetric replacement.
    identity(6 ** 6 * 5 + (n-2) * 5 ** 7 + 5 ** 6 * 6, 78125*n+170780)
    identity((5*n+1)*(78125*n+170780)-(625*n+671)**2, 91*(1025*n-3071))
    identity(91*(1025*(r+3)-3071), 93275*r+364)
    for dim in range(3, 101):
        lhs, rhs, cyc = sides([6] + [5] * (dim - 1), 6, 1)
        assert cyc == 78125*dim+170780
        assert lhs-rhs == 91*(1025*dim-3071) >= 364
    print('N(6,1)<=2: D_n(6,5,...,5)=91(1025n-3071)=93275(n-3)+364, checked n=3..100')
    lhs, rhs, cyc = sides([6, 5, 5], 6, 1)
    assert (cyc, lhs, rhs, lhs-rhs) == (405155, 6482480, 6482116, 364)
    print(f'COUNTEREXAMPLE n=3, x=(6,5,5), a=6,b=1,c=4: cyclic={cyc}, LHS={lhs}, RHS={rhs}, D={lhs-rhs}')

    # A second unbounded positive family, with positivity certified by coefficients.
    d3 = (3*m+1)*((m+1)**6*m+m**7+m**6*(m+1))-((m+1)**4+2*m**4)**2
    identity(d3, (3*m*m+3*m+1)*(m**3-4*m*m-4*m-1))
    identity((m**3-4*m*m-4*m-1).subs(m,r+5), r**3+11*r*r+31*r+4)
    for base in range(5, 101):
        lhs, rhs, _ = sides([base+1, base, base], 6, 1)
        assert lhs-rhs == (3*base*base+3*base+1)*(base**3-4*base*base-4*base-1) > 0
    print('FAMILY D_3(m+1,m,m)=(3m^2+3m+1)(m^3-4m^2-4m-1)>0 for integer m>=5')

    # Exact threshold values. The infinite n>=4 exclusion additionally uses
    # monotonicity of sin on [0,pi/2], explained in the README and draft.
    exact = [sp.simplify(48*sp.sin(sp.pi/sp.Integer(k))**2) for k in [2,3,4]]
    assert exact == [sp.Integer(48),sp.Integer(36),sp.Integer(24)]
    assert 36 <= exact[0] and 36 == exact[1] and exact[2] < 36
    print('N*(6,1)=3 threshold endpoints: n=2 -> 48, n=3 -> 36, n=4 -> 24, against 36')
    print('n>=4 analytic obligation: sin(pi/n)^2<=sin(pi/4)^2=1/2, hence 48 sin(pi/n)^2<=24<36')

    # Constant vector has a vanishing quadratic term but a positive cubic term.
    normalized = (3+t)*((1+t)**6+1+(1+t))-((1+t)**4+2)**2
    identity(normalized, -t**3*(t*t+3*t+3)*(t**3+4*t*t+4*t-1))
    assert sp.expand(normalized).coeff(t,2) == 0
    assert sp.expand(normalized).coeff(t,3) == 3
    assert normalized.subs(t,sp.Rational(1,5)) == sp.Rational(364,390625)
    print('MECHANISM D_3(1+t,1,1)=-t^3(t^2+3t+3)(t^3+4t^2+4t-1); quadratic=0, cubic=3')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
