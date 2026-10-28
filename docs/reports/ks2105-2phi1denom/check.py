#!/usr/bin/env python3
"""Exact finite checks of Kim--Stanton arXiv:2105.10096v3, Conjectures 8.3--8.4.

Python 3 and SymPy. No numeric approximation is used. These are finite checks,
not a proof for all n. Source fidelity: the Heine contiguous relation preceding
the lemma labelled lem:Heine, and both published continued fractions (labels lem:Heine,
prop:q-Norlund), evaluated independently of the coefficient quotient recurrence.
"""
from fractions import Fraction
import sys
import sympy
from sympy.polys.domains import QQ
from sympy.polys.fields import field


def poch(u, q, n):
    ans = q**0
    for j in range(n):
        ans *= 1 - u*q**j
    return ans


def series(a, b, c, q, degree):
    return [poch(a, q, n)*poch(b, q, n)/(poch(c, q, n)*poch(q, q, n))
            for n in range(degree+1)]


def quotient(g, f):
    assert f[0] == 1
    h = []
    for n in range(len(f)):
        h.append(g[n] - sum((f[i]*h[n-i] for i in range(1, n+1)), f[0]*0))
    return h


def mul(f, g):
    return [sum((f[i]*g[n-i] for i in range(n+1)), f[0]*0)
            for n in range(len(f))]


def reciprocal(f):
    return quotient([f[0]] + [f[0]*0]*(len(f)-1), f)


def heine(a, b, c, q, n):
    tail = [q**0] + [q*0]*n
    for m in range(n+2, 0, -1):
        if m % 2:
            r = (m-1)//2
            beta = (1-b*q**r)*(a-c*q**r)*q**r/((1-c*q**(2*r))*(1-c*q**(2*r+1)))
        else:
            r = m//2
            beta = (1-a*q**r)*(b-c*q**r)*q**(r-1)/((1-c*q**(2*r-1))*(1-c*q**(2*r)))
        tail = reciprocal([q**0] + [-beta*tail[j-1] for j in range(1, n+1)])
    return tail


def norlund(a, b, c, q, n):
    tail = [q*0]*(n+1)
    for m in range(n+1, -1, -1):
        bm = (a+b-a*b*q**m-a*b*q**(m+1))*q**m/(1-c*q**m)
        am = -(1-a*q**(m+1))*(1-b*q**(m+1))*c*q**m/((1-c*q**m)*(1-c*q**(m+1)))
        lm = (1-a*q**(m+1))*(1-b*q**(m+1))*a*b*q**(2*m+1)/((1-c*q**m)*(1-c*q**(m+1)))
        den = [q**0] + [-am*tail[j-1] - (lm*tail[j-2] if j >= 2 else 0)
                              - (bm if j == 1 else 0) for j in range(1, n+1)]
        tail = reciprocal(den)
    return tail


def denominator(c, q, n):
    ans = q**0
    for k in range(n+1):
        ans *= (1-c*q**k)**((n+1)//(k+1))
    return ans


def main():
    print('Python', sys.version.split()[0], 'SymPy', sympy.__version__, flush=True)
    K, a, b, c, q = field('a,b,c,q', QQ)
    nmax = 4
    f = series(a, b, c, q, nmax)
    g_single = series(a*q, b, c*q, q, nmax)
    shifted = series(a*q, b*q, c*q*q, q, nmax)
    beta = (1-b)*(a-c)/((1-c)*(1-c*q))
    for n in range(nmax+1):
        assert g_single[n]-f[n] == (beta*shifted[n-1] if n else 0)
    print('FIDELITY: published Heine contiguous relation, symbolic coefficients n=0..4: OK', flush=True)
    polynomials = {}
    for label, g in [('8.3', series(a*q, b*q, c*q, q, nmax)), ('8.4', g_single)]:
        gamma = quotient(g, f)
        counts = []
        polynomials[label] = []
        for n in range(nmax+1):
            h = denominator(c, q, n)*gamma[n]/(1-c)
            assert h.denom == K.ring.one, (label, n, h.denom)
            assert all(QQ.denom(coef) == 1 for coef in h.numer.values()), (label, n)
            polynomials[label].append(h)
            counts.append(len(h.numer))
            print('INTEGRAL:', label, 'n='+str(n), 'denominator=1', 'integer_coefficients=True',
                  'terms='+str(counts[-1]), flush=True)
        assert counts == ([1, 6, 34, 209, 678] if label == '8.3' else [1, 4, 24, 162, 530])
    grid = [(Fraction(2), Fraction(3), Fraction(5), Fraction(2)),
            (Fraction(-2), Fraction(3, 2), Fraction(-1, 3), Fraction(2, 3)),
            (Fraction(0), Fraction(0), Fraction(2, 5), Fraction(3, 2)),
            (Fraction(1), Fraction(7), Fraction(-3), Fraction(-2)),
            (Fraction(3, 4), Fraction(-2, 7), Fraction(5, 3), Fraction(1, 2))]
    for idx, (aa, bb, cc, qq) in enumerate(grid):
        degree = 12
        ff = series(aa, bb, cc, qq, degree)
        for label, gg, cf in [
            ('8.3', series(aa*qq, bb*qq, cc*qq, qq, degree), norlund),
            ('8.4', series(aa*qq, bb, cc*qq, qq, degree), heine)]:
            gamma = quotient(gg, ff)
            assert gamma == cf(aa, bb, cc, qq, degree), (idx, label)
            for n in range(nmax+1):
                symbolic = polynomials[label][n].as_expr()
                value = symbolic.subs(dict(zip(sympy.symbols('a b c q'),
                    [sympy.Rational(v.numerator, v.denominator) for v in (aa, bb, cc, qq)])))
                assert value == sympy.Rational(denominator(cc, qq, n)*gamma[n]/(1-cc))
        print('FIDELITY: exact rational grid', idx, 'both published continued fractions n=0..12: OK', flush=True)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
