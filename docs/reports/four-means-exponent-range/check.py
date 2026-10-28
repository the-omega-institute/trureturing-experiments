#!/usr/bin/env python3
"""Exact checks for Meštrović–Andjić Conjecture 1.6; no floating-point evidence."""
from fractions import Fraction as F
from math import comb
import sympy as S


def require(test, message):
    if not test:
        raise AssertionError(message)


def zero(expr):
    require(S.simplify(expr) == 0, str(expr))


def power_bounds(base, exponent, steps=112):
    """Enclose a positive rational base to a rational power by rational bisection."""
    base, exponent = F(base), F(exponent)
    require(base > 0, "positive base")
    if exponent < 0:
        lo, hi = power_bounds(base, -exponent, steps)
        require(lo > 0, "reciprocal bound")
        return 1 / hi, 1 / lo
    p, q = exponent.numerator, exponent.denominator
    target = base ** p
    if q == 1:
        return target, target
    lo, hi = F(0), max(F(1), target)
    for _ in range(steps):
        mid = (lo + hi) / 2
        if mid ** q <= target:
            lo = mid
        else:
            hi = mid
    require(lo ** q <= target <= hi ** q, "root enclosure certificate")
    return lo, hi


def difference_bounds(n, u):
    """F_n(u) = (2-u^2)^(n/2) + u^(2n) - 1 - u^n."""
    qlo, qhi = power_bounds(2 - u * u, n / 2)
    hlo, hhi = power_bounds(u, 2 * n)
    glo, ghi = power_bounds(u, n)
    return qlo + hlo - 1 - ghi, qhi + hhi - 1 - glo


# Fidelity: reproduce Theorem 1.1, (1) and the n=2 instance of (2).
x, y = S.symbols('x y', positive=True)
A, G2, H, Q2 = (x + y) / 2, x * y, 2 * x * y / (x + y), (x*x + y*y) / 2
zero(A*A*G2 - Q2*H*H - x*y*(x-y)**4/(4*(x+y)**2))
zero(Q2 + H*H - A*A - G2 - (x-y)**4/(4*(x+y)**2))
values = sorted({F(a, b) for a in range(1, 9) for b in range(1, 5)})
count = 0
for a in values:
    for b in values:
        ar, ge2, ha, qu2 = (a+b)/2, a*b, 2*a*b/(a+b), (a*a+b*b)/2
        require(ar*ar*ge2 >= qu2*ha*ha, "source (1), squared positive sides")
        require(qu2+ha*ha >= ar*ar+ge2, "source (2), n=2")
        require((qu2+ha*ha == ar*ar+ge2) == (a == b), "n=2 equality")
        count += 1
print(f'FIDELITY: source Theorem 1.1 (1), (2) at n=2: identities and {count} rational pairs OK')
print('SOURCE CAVEAT: at n=0 all four positive means have power 1; its stated equality iff x=y needs n!=0')

# Independently check the normalized identities and the local obstruction.
u, n = S.symbols('u n', real=True)
f = (2-u*u)**(n/2) + u**(2*n) - 1 - u**n
zero(f.subs(u, 1))
zero(S.diff(f, u).subs(u, 1))
zero(S.diff(f, u, 2).subs(u, 1) - 4*n*(n-1))
zero(2-u*u-(1+u-u*u)**2-(1-u)**3*(1+u))
zero((2-u*u)*u*u - (1-(1-u*u)**2))
t = S.symbols('t', real=True)
zero(1/t+t*t-1-t-(t-1)**2*(t+1)/t)
e = S.symbols('e', real=True)
eps_f = (1+e*e)**(n/2)+(1-e*e)**n-1-(1-e*e)**(n/2)
for order, expected in [(0, 0), (2, 0), (4, n*(n-1)/2)]:
    zero(S.diff(eps_f, e, order).subs(e, 0)/S.factorial(order)-expected)
print('LOCAL: F(1)=F\'(1)=0; F\'\'(1)=4*n*(n-1); epsilon^4 coefficient=n*(n-1)/2')
print('GLOBAL IDENTITIES: quadratic comparison, qu<=1 reduction, negative-power factorization OK')

# Exact infinite subfamily at n=1/2. For 1/2<=t<1, set u=t^2,
# x=1+sqrt(1-t^4), y=1-sqrt(1-t^4). Both are positive and unequal.
# Compare fourth powers of the positive quantities (2-t^4)^(1/4)
# and R=1+t-t^2. R^4-(2-t^4)=(t-1)^2*(t+1)*P(t).
P = t**5-3*t**4+4*t*t+3*t-1
zero((1+t-t*t)**4-(2-t**4)-(t-1)**2*(t+1)*P)
s = S.symbols('s', real=True)
polynomial = S.Poly(S.expand(P.subs(t, (1+s)/2)), s)
coeff = [polynomial.nth(j) for j in range(6)]
bernstein = [S.simplify(sum(coeff[j]*S.Rational(comb(k,j),comb(5,j))
                              for j in range(k+1))) for k in range(6)]
zero(polynomial.as_expr()-sum(bernstein[k]*comb(5,k)*s**k*(1-s)**(5-k)
                             for k in range(6)))
require(all(c.is_positive is True for c in bernstein), 'positive Bernstein coefficients')
print('HALF-POWER FAMILY: P((1+s)/2) Bernstein coefficients = '+str(bernstein))
print('HALF-POWER FAMILY: exact fourth-power gap positive for 1/2<=t<1; hence strict failures')

# The given algebraic anchor, with exact fourth-power comparison.
xa, ya = 2+S.sqrt(3), 2-S.sqrt(3)
zero(xa*ya-1)
zero((xa+ya)/2-2)
zero(2*xa*ya/(xa+ya)-S.Rational(1,2))
zero((xa*xa+ya*ya)/2-7)
zero((1+1/S.sqrt(2))**4-(S.Rational(17,4)+3*S.sqrt(2)))
require(F(18) > F(121,16), '3*sqrt(2)>11/4, positive sides squared')
print('ANCHOR: n=1/2, x=2+sqrt(3), y=2-sqrt(3): A=2,G=1,H=1/2,Q=sqrt(7); strict failure OK')

# A rational exponent grid checks both strict directions at actual positive pairs.
exponents = sorted({F(p,q) for q in range(2,13) for p in range(1,q)})
failures, reverse_failures = [], []
for exponent in exponents:
    for k in range(1,9):
        near = 1-F(1,10**k)
        lo, hi = difference_bounds(exponent, near)
        if hi < 0:
            failures.append((exponent,near,hi));break
    else:
        raise AssertionError(f'no exact negative enclosure for n={exponent}')
    for k in range(1,65):
        small = F(1,10**k)
        lo, hi = difference_bounds(exponent,small)
        if lo > 0:
            reverse_failures.append((exponent,small,lo));break
    else:
        raise AssertionError(f'no exact positive enclosure for n={exponent}')
print(f'EXACT GRID: {len(exponents)} rational exponents in (0,1), each with certified F<0 and F>0')
for exponent in [F(1,4),F(1,2),F(3,4)]:
    a=next(r for r in failures if r[0]==exponent)
    b=next(r for r in reverse_failures if r[0]==exponent)
    print(f'GRID WITNESS: n={exponent}, F({a[1]})<0 and F({b[1]})>0 (112 bisections, exact rationals)')
print('SCOPE: finite grid and symbolic identities support the prose argument; no Lean proof or novelty certification')
print('ALL_CHECKS_OK')
