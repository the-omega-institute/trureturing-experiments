#!/usr/bin/env python3
"""Exact checks for the six-letter three-ball maximum; no Lean proof."""
from collections import defaultdict
from itertools import product
from math import comb, factorial
import sympy as sp

x, y, z, a = sp.symbols('x y z a')
G = 1 + 5*x*y*z
F = (x+y*z+4*x*y*z, y+x*z+4*x*y*z, z+x*y+4*x*y*z)
P = sp.prod(F)
H = x*y+x*z+y*z+3*x*y*z
D = (1-x)*(1-y)*(1-z)
Q = 1-x-y-z-5*x*y-5*x*z-5*y*z-19*x*y*z
assert sp.expand(P-H**2*G-x*y*z*D*Q) == 0
print('Exact polynomial identity: P - H^2 G = xyz D Q')

# Fixed Laurent expansion: total surplus r+2u <= 4 is sufficient.
# This symbolic calculation is independent of the bounded integer enumeration below.
laurent = {(0, 0, 0): sp.Integer(1)}
for axis in range(3):
    out = defaultdict(lambda: sp.Integer(0))
    for r in range(5):
        for u in range((4-r)//2+1):
            e = [r+u]*3
            e[axis] = -r
            w = sp.prod(a-j for j in range(r+u))*4**u/sp.Integer(factorial(r)*factorial(u))
            for old, c in laurent.items():
                new = tuple(old[i]+e[i] for i in range(3))
                if sum(new) <= 4:
                    out[new] += c*w
    laurent = dict(out)
hq = sp.Poly(sp.expand(H*Q), x, y, z)
symbolic = sp.expand(sum(c*laurent.get(tuple(2-i for i in e), 0)
                         for e, c in hq.terms()))
assert sp.expand(symbolic - (12*a**3+9*a**2-156*a-57)) == 0
print('Universal symbolic coefficient in a=k-1:', symbolic)
h = sp.symbols('h')
cubic = lambda k: 12*k**3-27*k**2-138*k+96
assert sp.expand(cubic(5+h)-(12*h**3+153*h**2+492*h+231)) == 0
print('Positive-sign expansion at k=5+h: 12*h^3 + 153*h^2 + 492*h + 231')

# Sparse coordinate-by-coordinate enumerator, using only integer arithmetic.
def mul_bounded(poly, factor, bound):
    out = defaultdict(int)
    for e, c in poly.items():
        for f, w in factor.items():
            g = tuple(e[i]+f[i] for i in range(3))
            if max(g) <= bound:
                out[g] += c*w
    return {e:c for e,c in out.items() if c}

def factor_dict(expr):
    return {e:int(c) for e,c in sp.Poly(sp.expand(expr), x, y, z).terms()}

def core(q, aa, b, bound, common=0):
    fs = [
        {(1,0,0):1, (0,1,1):1, (1,1,1):q-2},
        {(0,1,0):1, (1,0,1):1, (1,1,1):q-2},
        {(0,0,1):1, (1,1,0):1, (1,1,1):q-2}]
    hh = {(1,1,0):1, (1,0,1):1, (0,1,1):1, (1,1,1):q-3}
    gg = {(0,0,0):1, (1,1,1):q-1}
    poly = {(0,0,0):1}
    for f, count in [(f,aa) for f in fs]+[(hh,b),(gg,common)]:
        for _ in range(count):poly = mul_bounded(poly, f, bound)
    return poly

def ball_count(poly, t):
    return sum(c for e,c in poly.items() if max(e) <= t)

coefficients = []
for k in range(1,16):
    bound = k+2
    pa = core(6,k,1,bound)
    pb = core(6,k-1,3,bound,common=1)
    assert ball_count(pa,k+1) == ball_count(pb,k+1) == 6
    assert all(ball_count(pa,r)-ball_count(pb,r)==0 for r in range(k+2))
    cpoly = core(6,k-1,1,bound)
    cpoly = mul_bounded(cpoly, factor_dict(x*y*z*Q), bound)
    assert all(cpoly.get((r,r,r),0)==0 for r in range(k+2))
    got = cpoly.get((bound,bound,bound),0)
    assert got == cubic(k) == ball_count(pa,bound)-ball_count(pb,bound)
    # Independent convolution check with nonconstant N and a larger radius.
    t = k+3
    cc = core(6,k-1,1,t)
    cc = mul_bounded(cc,factor_dict(x*y*z*Q),t)
    for n_extra in (0,1,3):
        lhs = ball_count(core(6,k,1,t,common=n_extra),t)-ball_count(core(6,k-1,3,t,common=n_extra+1),t)
        rhs = sum(comb(n_extra,j)*5**j*cc.get((t-j,)*3,0)
                  for j in range(min(n_extra,t)+1))
        assert lhs == rhs
    coefficients.append(got)
    print(f'k={k}: minimal tie=6; all c_r=0 through r={k+1}; c_(k+2)={got}; convolution OK')
assert coefficients[:5] == [-57,-192,-237,-120,231]

# Exhaustive enumeration of actual words, independent of polynomial counting.
def centers(aa,b,n):
    columns = [(1,0,0)]*aa+[(0,1,0)]*aa+[(0,0,1)]*aa+[(0,1,2)]*b+[(0,0,0)]*(n-3*aa-b)
    assert len(columns)==n
    return tuple(tuple(col[i] for col in columns) for i in range(3))

def count_words(q, cs, t):
    return sum(all(sum(u!=v for u,v in zip(w,c))<=t for c in cs)
               for w in product(range(q),repeat=len(cs[0])))
for q in (3,4,5,6,7):
    ca,cb = centers(1,1,4),centers(0,3,4)
    assert all(sum(u!=v for u,v in zip(c,d))==3 for c,d in ((ca[0],ca[1]),(ca[0],ca[2]),(ca[1],ca[2])))
    got_a,got_b = count_words(q,ca,2),count_words(q,cb,2)
    # Source v2, proof of Theorem 25, equations |C|d and |C|d+1.
    assert (got_a,got_b)==(q,6)
    print(f'Fidelity q={q}, d=3, n=4, radius 2: A={got_a}, B={got_b} (source q and 6)')
ca,cb = centers(1,1,4),centers(0,3,4)
assert (count_words(6,ca,3),count_words(6,cb,3))==(189,246)
print('Exhaustive 6^4 words, d=3, t=3: A=189, B=246, A-B=-57')
print('ALL_CHECKS_OK')
