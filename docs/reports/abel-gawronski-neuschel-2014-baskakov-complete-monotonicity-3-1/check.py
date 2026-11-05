#!/usr/bin/env python3
"""Exact certificates for Abel--Gawronski--Neuschel, arXiv:1411.7945, Conjecture 3.1.

Run with Python 3 and SymPy. No approximate numbers participate in any check.
The continuous-parameter interval and eventual derivative signs require the
analytic arguments in the preregistration; this program checks their finite inputs.
"""
from fractions import Fraction
from math import comb
import sympy as s

x, y, a, u = s.symbols('x y a u', real=True)
half = s.Rational(1, 2)
sqrt2 = s.sqrt(2)


def zero(expr):
    assert s.cancel(s.expand(expr), extension=True) == 0


def reciprocal_coefficients(d, order):
    coeff = [Fraction(1, d[0])]
    for m in range(1, order + 1):
        coeff.append(-sum((d[j] * coeff[m-j]
                          for j in range(1, min(m, len(d)-1)+1)), Fraction()) / d[0])
    return coeff


# Source Theorem 2.3 (eq1), n=1,2,3: independently sum the squared
# Baskakov series with its generating polynomial in k, then compare the
# source's finite positive-coefficient formula preceding Remark 2.4.
z = s.symbols('z')
k = s.symbols('k', integer=True, nonnegative=True)
for n in range(1, 4):
    pk = s.prod(k+j for j in range(1, n)) / s.factorial(n-1)
    polynomial = s.Poly(s.expand(pk**2), k)
    moments = [1/(1-z)]
    for j in range(1, polynomial.degree()+1):
        moments.append(s.cancel(z*s.diff(moments[-1], z)))
    generating = sum(polynomial.nth(j)*moments[j]
                     for j in range(polynomial.degree()+1))
    direct = generating.subs(z, (x/(1+x))**2)/(1+x)**(2*n)
    source = (s.Rational(comb(2*n-2,n-1), 2**(2*n-2)) *
              sum(s.Rational(comb(n-1,j)**2, comb(2*n-2,2*j)) /
                  (1+2*x)**(2*j+1) for j in range(n)))
    zero(direct-source)
    print(f'FIDELITY source squared-Baskakov finite identity n={n}: exact equality')

# alpha=1: the literal binomial series is geometric for every even r.
for r in (2, 4, 6, 8):
    for j in range(20):
        assert s.binomial(-1, j)**r == 1
    zero((1+x)**(-r)/(1-(x/(1+x))**r) - 1/((1+x)**r-x**r))
    print(f'GEOMETRIC alpha=1 r={r}: exact rational identity')

# Rational Laplace transforms of the damped cosine kernels, a=x+1/2.
transforms = {
    2: [(half, s.Integer(0))],
    4: [(s.Integer(1), s.Integer(0)), (-s.Integer(1), half)],
    6: [(s.Rational(8,3), s.Integer(0)),
        (-s.Integer(3), 1/(2*s.sqrt(3))), (s.Rational(1,3), s.sqrt(3)/2)],
    8: [(s.Integer(8), s.Integer(0)), (s.Integer(2), half),
        (-(5+7*sqrt2/2), (sqrt2-1)/2),
        (-(5-7*sqrt2/2), (sqrt2+1)/2)]
}
for r, terms in transforms.items():
    laplace = sum(A*a/(a*a+b*b) for A,b in terms)
    zero(laplace.subs(a,x+half)-1/((1+x)**r-x**r))
    print(f'KERNEL r={r}: exact rational Laplace identity')
# Nonnegative factorization for r=6 after cos(3u)=4cos(u)^3-3cos(u).
zero(s.Rational(8,3)-3*u+s.Rational(1,3)*(4*u**3-3*u)
     -s.Rational(4,3)*(u-1)**2*(u+2))
print('SURVIVING r=2,4,6: nonnegative densities via 1/2, 1-cos(u), (4/3)(cos(u)-1)^2(cos(u)+2)')

# Exact algebraic inputs for the strict trigonometric sign certificate.
# 7 < 5sqrt(2) < 29/4 => 0 < (5sqrt(2)-7)pi < pi/4.
assert 7**2 < 5**2*2 < s.Rational(29,4)**2
assert 6**2 < 5**2*2
A = 5+7*sqrt2/2
B = 5-7*sqrt2/2
zero(8-2+A*u+B*u-(6+10*u))
print('NEGATIVE KERNEL t=10*pi: bracket = 6+10*cos(5*sqrt(2)*pi) < 6-5*sqrt(2) < 0')

# Rational Taylor jet of F(2+y), two independent exact recurrences.
q = s.Poly(s.expand((3+y)**8-(2+y)**8),y)
d = [int(q.nth(j)) for j in range(8)]
assert d == [6305,16472,18620,11816,4550,1064,140,8]
c = reciprocal_coefficients(d,75)
# Independent Leibniz recursion: q F = 1 in ordinary derivative units.
jet = [Fraction(1,d[0])]
for m in range(1,76):
    jet.append(-sum((comb(m,j)*s.factorial(j)*d[j]*jet[m-j]
                    for j in range(1,min(m,7)+1)), Fraction())/d[0])
for m in range(76):
    assert Fraction(jet[m]) == s.factorial(m)*c[m]
    assert sum((d[j]*c[m-j] for j in range(min(m,7)+1)), Fraction()) == (1 if m==0 else 0)
assert c[75] > 0
assert c[75].numerator == 52526631898460004553579116030127113870398461748869136952441542820718397538730226967476064508798375739888524672284671945933858519867009460198016371252290979899412571544360280998730562833801497863350200902760105459773399085006189628378378200422549498764460032
print('JET denominator coefficients:', ','.join(map(str,d)))
print('JET (-1)^75 F^(75)(2)/75! =', -c[75])
print('JET exact sign: negative (two recurrences agree through order 75)')

# Nondegenerate infinite Baskakov parameter family: for every positive n,
# c=n gives n/c=1 and psi_{n,n}^[8](x)=F(n*x).
# Symbolic identity in n proves that the finite samples do not define the scope.
scale = s.symbols('scale', positive=True)
zero((1+scale*x)**(-8)/(1-(scale*x/(1+scale*x))**8)
     -1/((1+scale*x)**8-(scale*x)**8))
for n in range(1,9):
    dn = [d[j]*n**j for j in range(8)]
    cn = reciprocal_coefficients(dn,75)
    assert cn[75] == n**75*c[75] > 0
    print(f'FAMILY n=c={n}, x=2/{n}: positive derivative 75 (exact)')
print('FAMILY symbolic n=c>0 identity checked; integer n>=1 yields infinitely many counterexamples')
print('INTERVAL alpha near 1: analytic uniform-differentiation obligation; no numeric epsilon claimed')
print('ALL_CHECKS_OK')
