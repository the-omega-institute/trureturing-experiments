#!/usr/bin/env python3
"""Exact bounded checks for Voloshyn, arXiv:2609.31189v1, Conjecture 8.1.

Requires Python 3 and SymPy. These checks are neither an all-order proof nor
a certification that a computed norm is a minimal polynomial.
"""
from fractions import Fraction as Q
import sympy as sp

x, y, z, t, w = sp.symbols('x y z t w')


def zero(expr):
    assert sp.cancel(expr) == 0, sp.factor(expr)


def mul(a, b, n):
    return [sum((a[j]*b[i-j] for j in range(i+1)), Q(0)) for i in range(n+1)]


def power(a, exponent, n):
    out = [Q(1)] + [Q(0)]*n
    for _ in range(exponent):
        out = mul(out, a, n)
    return out


def inverse(a, n):
    assert a[0] != 0
    out = [1/a[0]] + [Q(0)]*n
    for i in range(1, n+1):
        out[i] = -sum((a[j]*out[i-j] for j in range(1, i+1)), Q(0))/a[0]
    return out


def implicit_branch(q, yy, n):
    """Solve W(z,s)=0 at z(0)=y directly in truncated Q[[s]]."""
    zz = [yy] + [Q(0)]*n
    slope = (1-yy)**q * q * yy**(q-1)
    assert slope != 0
    for k in range(1, n+1):
        one_minus = [1-zz[0]] + [-a for a in zz[1:]]
        zq = power(zz, q, n)
        zq[0] -= yy**q
        residual = mul(power(one_minus, q, n), zq, n)[k]
        residual -= power(zz, q, n)[k-1]
        zz[k] = -residual/slope
    one_minus = [1-zz[0]] + [-a for a in zz[1:]]
    zq = power(zz, q, n)
    zq[0] -= yy**q
    residual = mul(power(one_minus, q, n), zq, n)
    for k in range(1, n+1):
        residual[k] -= power(zz, q, n)[k-1]
    assert all(a == 0 for a in residual)
    return zz


def branch_f(q, yy, n):
    zz = implicit_branch(q, yy, n)
    one_minus = [1-zz[0]] + [-a for a in zz[1:]]
    invz = power(inverse(zz, n), q-1, n)
    bracket = [(q-1)*a-q*yy**q*b for a, b in zip(zz, invz)]
    bracket[0] += 1
    vv = mul(power(one_minus, q-1, n), bracket, n)
    vv[0] -= (1-yy)**q
    return [Q((-1)**(q+1), q**3)*a for a in vv]


def boundary_norm(q, yy):
    """Compute the norm in Q(y)[Z]/((1-Z)^q (Z^q-y^q))."""
    c = (1-yy)**q
    W = (-1)**q*(1-z)**q*(z**q-yy**q)
    fnum = sp.Rational((-1)**(q+1), q**3)*(
        (1-z)**(q-1)*(z**(q-1)*(1+(q-1)*z)-q*yy**q)-c*z**(q-1))
    numerator = sp.resultant(W, t*z**(q-1)-fnum, z)
    denominator = sp.resultant(W, z**(q-1), z)
    assert denominator != 0
    return sp.Poly(sp.cancel(numerator/denominator), t).as_expr()


def conjectured_boundary(q, xx):
    """Source (95)--(97), eliminating roots of unity by a resultant."""
    beta = sp.Rational((-1)**q, q**3)*(1+xx)**q
    beta_w = sp.Rational((-1)**q, q**3)*((1+xx)**q-(1+w*xx)**q)
    # zeta=1 contributes the source's separate t factor.
    residual_with_zero = sp.resultant(w**q-1, t-beta_w, w)
    return sp.expand((t-beta)**q*residual_with_zero)


def main():
    print('Arithmetic: integers, fractions.Fraction, SymPy exact rationals/polynomials')
    print('SymPy version:', sp.__version__)
    for q in range(2, 7):
        c = sp.Rational(q-1, q)
        numerator, exponent = c*x, 1
        rows = [numerator/(1+x)**exponent]
        for k in range(1, 4):
            lam = sp.Rational((q*k+q)*sp.prod(q*k+1-j for j in range(q)), q*q)
            # Differentiate N/(1+x)^e exactly, keeping the denominator factored.
            for _ in range(q+1):
                numerator = sp.expand(sp.diff(numerator, x)*(1+x)-exponent*numerator)
                exponent += 1
            numerator = sp.expand(x*numerator/lam)
            rows.append(numerator/(1+x)**exponent)
        zero(rows[1]-(-1)**q*sp.Rational(q-1, 2)*x/(1+x)**(q+2))
        zero(rows[2]-sp.Rational((q-1)*q*(q+1), 6)*x*(1-x)/(1+x)**(2*q+3))
        print(f'FIDELITY q={q}: source R1 and R2 reproduced', flush=True)
        for k, row in enumerate(rows):
            expanded = sp.Poly(sp.series(row, x, 0, 9).removeO(), x)
            for n in range(1, 9):
                predicted = (c*(-1)**(n+q*k-1)*q**(2*k)*sp.rf(sp.Rational(n, q), k)
                             * sp.factorial(n+q*k)
                             / (sp.factorial(n)*sp.factorial(k+1)*sp.factorial(q*k+1)))
                zero(expanded.nth(n)-predicted)
        # Solve the branch directly, independently of the Lagrange formula.
        for xx in [Q(1), Q(2), Q(3)]:
            coeffs = branch_f(q, -xx, 4)
            assert coeffs[0] == 0
            for k, row in enumerate(rows):
                expected = sp.cancel(row.subs(x, sp.Rational(xx.numerator, xx.denominator)))
                actual = coeffs[k+1]*((-1)**q*q*q)**(k+1)
                zero(sp.Rational(actual.numerator, actual.denominator)-expected)
        U = (1-z)**q*(1-y**q/z**q)
        V = (1-z)**(q-1)*(1+(q-1)*z-q*y**q/z**(q-1))
        zero(sp.diff(V, z)-(q-1)*z*sp.diff(U, z)/(1-z))
        zero(V.subs(z, y)-(1-y)**q)
        zero(V.subs(z, 1))
        num = sp.cancel(V.subs(z, y*w)-(1-y*w)**q).as_numer_denom()[0]
        zero(sp.rem(num, w**q-1, w))
        for expr, degree in [(U, 2*q), (V, 2*q-1)]:
            num, den = sp.fraction(sp.cancel(expr))
            pp = sp.Poly(num, z, domain=sp.QQ.frac_field(y))
            dd = sp.Poly(den, z, domain=sp.QQ.frac_field(y))
            assert sp.gcd(pp, dd).degree() == 0
            assert max(pp.degree(), dd.degree()) == degree
        assert sp.gcd(2*q, 2*q-1) == 1
        norm = boundary_norm(q, -x)
        zero(norm-conjectured_boundary(q, x))
        assert sp.Poly(norm, t).LC() == 1
        assert sp.degree(norm, t) == 2*q
        print(f'IDENTITY q={q}: coefficients, branch, primitive, degrees, symbolic boundary norm OK', flush=True)
    p4 = t*(27*t+(1+x)**3)**3*(27*t**2+9*x*(1+x)*t+x*x*(x*x+x+1))/3**12
    zero(boundary_norm(3, -x)-p4)
    p6_at_one = t*(125*t+32)**5*(125*t*t+75*t+11)**2/5**21
    zero(boundary_norm(5, -sp.Integer(1))-p6_at_one)
    zero(conjectured_boundary(1, x)-t*(t+1+x))
    print('FIDELITY: p=4 boundary, p=6 at x=1, and separately normalized Narayana endpoint OK')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
