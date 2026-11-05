#!/usr/bin/env python3
"""Exact preregistration checks for Karp, arXiv:1203.1482v1, Conjecture 2.

No floating arithmetic or numerical root finder is used. Routh's rule is
applied only when every pivot is nonzero and no all-zero row occurs.
The finite checks do not establish the quantified asymptotic family.
"""
from fractions import Fraction
from functools import reduce
from math import comb, gcd
import sys
import sympy as sp


def mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            c[i + j] += ai * bj
    return c


def poch(shift, n):
    out = [1]
    for j in range(n):
        out = mul(out, [shift + j, 1])
    return out


def q_integer(m, A):
    """A^(2m^2) Q_(2m)^(A,A), f_k=A^(-(k-m)^2), in Z[x]."""
    n = 2 * m
    ps = {c: [poch(c * A, k) for k in range(n + 1)] for c in (0, 1, 2)}
    out = [0] * (n + 1)
    for k in range(n + 1):
        w = comb(n, k) * A ** (2 * m * m - 2 * (k - m) ** 2)
        a = mul(ps[1][k], ps[1][n - k])
        b = mul(ps[2][k], ps[0][n - k])
        for j in range(n + 1):
            out[j] += w * (a[j] - b[j])
    assert out[-2:] == [0, 0], 'two leading coefficients must cancel'
    return out[:-2]


def primitive(row):
    d = reduce(gcd, (abs(x) for x in row), 0)
    assert d, 'all-zero Routh row: exceptional case is not implemented'
    return [x // d for x in row]


def routh(coeff):
    """Exact Routh sign sequence using only positive row rescalings."""
    a = list(reversed(coeff))
    width = (len(a) + 1) // 2
    rows = [primitive(a[::2] + [0] * (width - len(a[::2]))),
            primitive(a[1::2] + [0] * (width - len(a[1::2])))]
    for _ in range(len(a) - 2):
        p, q = rows[-2:]
        assert q[0] != 0, 'zero pivot: exceptional case is not implemented'
        # Usual new row = numerator / q[0]. Multiplication by |q[0]|
        # clears that denominator while preserving every sign.
        sg = 1 if q[0] > 0 else -1
        row = [sg * (q[0] * p[j + 1] - p[0] * q[j + 1])
               for j in range(width - 1)] + [0]
        rows.append(primitive(row))
    assert all(row[0] for row in rows)
    signs = [1 if row[0] > 0 else -1 for row in rows]
    changes = sum(a != b for a, b in zip(signs, signs[1:]))
    return ''.join('+' if s > 0 else '-' for s in signs), changes


def fidelity():
    x, alpha, beta = sp.symbols('x alpha beta')
    # Source Section 2, Lemma 3 (label lm:Chu): f_k=1 => Q_n identically 0.
    for n in range(9):
        q = sum(comb(n, k) * (sp.rf(x + alpha, k) * sp.rf(x + beta, n-k)
                - sp.rf(x + alpha + beta, k) * sp.rf(x, n-k)) for k in range(n+1))
        assert sp.Poly(sp.expand(q), x, alpha, beta).is_zero
    print('FIDELITY: source Lemma 3, Q_n=0 for f_k=1, symbolic alpha,beta,x; n=0..8')
    # Separate determinant/series convolution calculation (equations (2),(3)).
    t = sp.symbols('t')
    f = [1, 3, 7, 12, 14, 9]
    def series(shift):
        return sum(f[k] * sp.rf(x + shift, k) * t**k / sp.factorial(k) for k in range(6))
    det = sp.Poly(sp.expand(series(alpha)*series(beta)-series(alpha+beta)*series(0)), t)
    for n in range(6):
        direct = sum(comb(n,k)*f[k]*f[n-k]*(sp.rf(x+alpha,k)*sp.rf(x+beta,n-k)
                  - sp.rf(x+alpha+beta,k)*sp.rf(x,n-k)) for k in range(n+1))
        assert sp.expand(sp.factorial(n)*det.nth(n)-direct) == 0
    print('FIDELITY: determinant coefficient extraction equals source equation (3); n=0..5')


def limiting_checks():
    z, t = sp.symbols('z t')
    for m in range(2, 17):
        F = sp.Poly((z+1)**(2*m)-(z*(z+2))**m, z)
        assert F.degree() == 2*m-2 and F.LC() == m
        w = z*(z+2)
        assert sp.expand(F.as_expr()-sum(comb(m,j)*w**j for j in range(m))) == 0
        assert all(c > 0 for c in F.all_coeffs())
        # The central scaled summand at t=1/A tends exactly to F.
        c = sp.prod(z+1+j*t for j in range(m))**2
        d = sp.prod(z+2+j*t for j in range(m))*sp.prod(z+j*t for j in range(m))
        assert sp.expand((c-d).subs(t,0)-F.as_expr()) == 0
        if m >= 10:
            signs, count = routh([int(F.nth(j)) for j in range(2*m-1)])
            assert count >= 2*(m//10)
            print(f'LIMIT: m={m}, degree={F.degree()}, positive coefficients, Routh={signs}, RHP={count}')
    # Exact m=10 cyclotomic root location: Re(r)>0 iff sin(pi/10)<1/3.
    s = (sp.sqrt(5)-1)/4
    assert sp.simplify(sp.sin(sp.pi/10)-s) == 0
    assert s > 0 and s < sp.Rational(1,3)
    assert sp.sqrt(5) < sp.Rational(7,3)  # squared inequality 5 < 49/9
    assert 5 < sp.Rational(49,9)
    # Algebraic root identity, independently of choosing a square-root branch.
    u, v = sp.symbols('u v')
    F10 = (z+1)**20-(z*(z+2))**10
    replaced = sp.Poly(sp.expand(F10.subs(z,u-1)),u)
    even = sum(replaced.nth(2*j)*v**j for j in range(11))
    cyclo = sp.cyclotomic_poly(10,t)
    rational = sp.together(even.subs(v,1/(1-t)))
    numerator, denominator = sp.fraction(rational)
    assert sp.rem(sp.Poly(numerator,t),sp.Poly(cyclo,t)).is_zero
    assert sp.rem(sp.Poly(denominator,t),sp.Poly(cyclo,t)) != sp.Poly(0,t)
    print('LIMIT: m=10 cyclotomic root identity and strict positive real part checked exactly')


def family_instances():
    # Validate the sign-change implementation on known exact-root products.
    for coeff, expected in [([1,3,2],0),([-2,-1,1],1),([2,-2,1],2),([2,2,1],0)]:
        assert routh(coeff)[1] == expected
    print('ROUTH: four exact-root control polynomials checked')
    for m, A, expected in [(10,100,0), (10,1000,2), (10,10000,2), (11,1000,2), (12,1000,2)]:
        f = [Fraction(1,A**((k-m)**2)) for k in range(2*m+1)]
        assert all(x > 0 for x in f)
        assert all(f[k]**2/(f[k-1]*f[k+1]) == A**2 for k in range(1,2*m))
        q = q_integer(m,A)
        assert all(c > 0 for c in q)
        signs, count = routh(q)
        assert count == expected
        label = 'WITNESS' if count else 'FINITE_PARAMETER_CONTROL'
        print(f'{label}: m={m}, n={2*m}, alpha=beta={A}, degree={len(q)-1}, '
              f'positive coefficients, strict log-concavity ratio={A**2}, Routh={signs}, RHP={count}')
    print('SCOPE: finite witnesses and limiting identities; the universal large-A family is a prose proof obligation')


def main():
    print('Python',sys.version.split()[0], 'SymPy',sp.__version__, 'arithmetic=exact')
    fidelity()
    limiting_checks()
    family_instances()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
