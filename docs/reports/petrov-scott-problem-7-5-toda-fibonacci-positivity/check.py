#!/usr/bin/env python3
"""Exact algebraic checks for Petrov--Scott Problem 7.5; no numerical ODE evidence.

Fidelity: source (7.3) Motzkin moments and Example 7.3 Toda orientation.
Candidate: scaling, strict initial positivity, the finite B identity, coordinate
vector field, transverse derivative, zero eigenvector and atom-normalization
identities. Infinite spectral existence and uniform-in-word exit are hand
arguments, not consequences claimed from this finite computation.
"""
from fractions import Fraction as F
from itertools import product
from math import factorial
import sympy as S


def equal(a, b):
    d = S.cancel(S.expand(a - b))
    assert d == 0, (a, b, d)
    assert not S.sympify(a).has(S.Float)
    assert not S.sympify(b).has(S.Float)


def tri(diag, upper):
    n = len(diag)
    return S.Matrix(n, n, lambda i, j:
                    diag[i] if i == j else
                    upper[i] if j == i + 1 else
                    1 if i == j + 1 else 0)


def detq(diag, upper):
    """Rational continuant, independent of the symbolic matrix determinant."""
    d0, d1 = F(1), F(diag[0]) if diag else F(1)
    for i in range(1, len(diag)):
        d0, d1 = d1, diag[i] * d1 - upper[i-1] * d0
    return d1


def Aq(x, y, n):
    return detq([x(k) for k in range(1, n+1)],
                [y(k) for k in range(1, n)])


def Bq(x, y, k, r=0):
    return detq([y(r+1)] + [x(r+i+2) for i in range(1, k+1)],
                ([x(r+1)*y(r+2)] + [y(r+i+2) for i in range(1, k)]) if k else [])


def sq(x, y, w):
    k = 0
    while k < len(w) and w[k] == 1:
        k += 1
    if k == len(w):
        return Aq(x, y, k)
    assert w[k] == 2
    u = w[k+1:]
    return Bq(x, y, k, sum(u)) * sq(x, y, u)


def words(weight):
    if weight == 0:
        yield ()
    for a in (1, 2):
        if weight >= a:
            for w in words(weight-a):
                yield (a,) + w


# Reproduce the four normalized moment polynomials in source equation (7.3).
xs = S.symbols('x0:7')
ys = S.symbols('y0:7')
state = {0: S.Integer(1)}
moments = [S.Integer(1)]
for n in range(1, 5):
    nxt = {}
    for h, value in state.items():
        nxt[h] = nxt.get(h, 0) + value * xs[h+1]
        nxt[h+1] = nxt.get(h+1, 0) + value * ys[h+1]
        if h:
            nxt[h-1] = nxt.get(h-1, 0) + value
    state = nxt
    moments.append(S.expand(state.get(0, 0)))
source_moments = [1, xs[1], xs[1]**2+ys[1],
    xs[1]**3+2*xs[1]*ys[1]+xs[2]*ys[1],
    xs[1]**4+3*xs[1]**2*ys[1]+ys[1]**2+2*xs[1]*xs[2]*ys[1]
    +xs[2]**2*ys[1]+ys[1]*ys[2]]
for a, b in zip(moments, source_moments):
    equal(a, b)
print('FIDELITY: source (7.3), all four moment polynomials agree exactly')

# Example 7.3: q = exp(varrho), so time derivative is q*d/dq.
q = S.Symbol('q', positive=True)
for n in range(1, 9):
    xn, yn, yp = n+q-1, n*q, (n-1)*q
    equal(q*S.diff(xn, q), yn-yp)
    equal(q*S.diff(yn, q), yn*((n+1+q-1)-xn))
print('FIDELITY: Example 7.3 Toda equations, n=1..8, correct backward-time orientation')

# Strict positivity is checked for the original clone word values, not only
# non-strict normalized boundary inequalities.
u, t = lambda k: F(k), lambda k: F(k)
for n in range(13):
    assert Aq(u, t, n) == 1
for r in range(6):
    for k in range(8):
        assert Bq(u, t, k, r) == r+1
all_words = [w for n in range(13) for w in words(n)]
assert all(sq(u, t, w) > 0 for w in all_words)
print(f'INITIAL: A_n=1 (n=0..12); B_k(shift r)=r+1 (k=0..7,r=0..5); {len(all_words)} positive words (weight<=12)')

# All independent choices in a rational tail grid: c_1,...,c_4=1/k,
# c_k in {1/(2k),1/k,3/(2k)} for k=5..8.
small_words = [w for n in range(9) for w in words(n)]
checks = 0
for tail in product((F(1, 2), F(1), F(3, 2)), repeat=4):
    cs = {k: F(1,k) for k in range(1,5)}
    cs.update({k: tail[k-5]/k for k in range(5,9)})
    x = lambda k: k*cs[k]
    y = lambda k: k*cs[k]*cs[k+1]
    assert cs[1]-4*cs[2]+8*cs[3]-6*cs[4] == F(1,6)
    for w in small_words:
        scale = S.prod(cs[k] for k in range(1,sum(w)+1))
        assert sq(x,y,w) == scale*sq(u,t,w) > 0
        checks += 1
print(f'FAMILY: 81 independent rational tails, {checks} exact scaling/positivity comparisons (weight<=8), D=1/6')

# Generic finite clone determinant identity, with actual tridiagonal matrices.
ts = [S.Integer(0)] + list(S.symbols('t1:13'))
for r in range(4):
    for N in range(6):
        k = N+1
        diag = [ts[r+1]] + [1+ts[r+i+1] for i in range(1,k+1)]
        upper = [(1+ts[r])*ts[r+2]] + [ts[r+i+2] for i in range(1,k)]
        actual = tri(diag, upper).det(method='domain-ge')
        series = 1 + sum(S.prod(ts[i] for i in range(r+3,r+3+j)) for j in range(1,N+1))
        expected = ts[r+1] + (ts[r+1]-ts[r]-1)*ts[r+2]*series
        equal(actual, expected)
print('DETERMINANTS: generic B_(N+1)(shift r) identity, N=0..5,r=0..3 (24 symbolic matrix determinants)')

# Generic coordinate vector field: verify induced Toda equations before
# substituting the saturated boundary t_k=k.
cs = [S.Integer(0)] + list(S.symbols('c1:13', positive=True))
def cflow(n):
    return cs[n]*(cs[n+1]*ts[n]-cs[n]*ts[n-1])
def tflow(n):
    return ts[n]*(cs[n+1]*(1+ts[n])-cs[n]-cs[n+2]*ts[n+1]) if n else 0
for n in range(1,7):
    x = cs[n]*(1+ts[n-1])
    y = cs[n]*cs[n+1]*ts[n]
    xp = cflow(n)*(1+ts[n-1]) + cs[n]*tflow(n-1)
    yp = cflow(n)*cs[n+1]*ts[n]+cs[n]*cflow(n+1)*ts[n]+cs[n]*cs[n+1]*tflow(n)
    equal(xp, y-(cs[n-1]*cs[n]*ts[n-1] if n>1 else 0))
    equal(yp, y*(cs[n+1]*(1+ts[n])-x))
D = cs[1]-4*cs[2]+8*cs[3]-6*cs[4]
equal((tflow(2)-tflow(1)).subs({ts[i]:i for i in range(1,5)}),D)
L = S.Symbol('L', positive=True)
equal(D.subs({cs[i]:L/i for i in range(1,5)}), L/6)
print('VECTOR FIELD: generic induced Toda equations n=1..6; delta_prime=D; c_k=L/k gives D=L/6')

# P_n(0), the monic norm product and p_n(0)^2, with all c and t free.
p0, p1 = S.Integer(1), -cs[1]
for n in range(1,9):
    pn = p1
    equal(pn, (-1)**n*S.prod(cs[i] for i in range(1,n+1)))
    norm_sq = S.prod(cs[i]*cs[i+1]*ts[i] for i in range(1,n+1))
    Tn = S.prod(ts[i] for i in range(1,n+1))
    equal(pn**2/norm_sq, cs[1]/(cs[n+1]*Tn))
    p0, p1 = p1, -cs[n+1]*(1+ts[n])*p1 - cs[n]*cs[n+1]*ts[n]*p0
print('ATOM OBSTRUCTION: generic p_n(0)^2=c_1/(c_(n+1)*T_n), n=1..8 (essential c ratio retained)')

# The proposed zero eigenvector, checked in each row using exact radicals.
v = lambda k: (-1)**(k-1)*S.sqrt(cs[1]/(cs[k]*S.factorial(k-1)))
for k in range(1,9):
    row = k*cs[k]*v(k) + S.sqrt(k*cs[k]*cs[k+1])*v(k+1)
    if k>1:
        row += S.sqrt((k-1)*cs[k-1]*cs[k])*v(k-1)
    equal(S.simplify(row),0)
n = S.Symbol('n', integer=True, nonnegative=True)
equal(S.summation((n+1)/S.factorial(n),(n,0,S.oo)), 2*S.E)
print('ZERO ATOM: Jv=0 rows 1..8 for generic positive c; exact sum (n+1)/n! = 2*e, so explicit-family mass=1/(2*e)')
print('SCOPE: finite exact algebra only; infinite spectral existence, tilting and immediate exit require the hand argument and later Lean proof')
print('ALL_CHECKS_OK')
