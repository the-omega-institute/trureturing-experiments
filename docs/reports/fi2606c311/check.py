#!/usr/bin/env python3
"""Exact algebra for Frank--Ivkovic arXiv:2606.16597v3, Conjecture 3.11.

This checks identities and rational samples, not the analytic bidual bridge or a Lean proof.
Requires Python >= 3.12 and SymPy. No floating-point arithmetic is used.
"""
import sys
from fractions import Fraction as Q
from itertools import product
import sympy as s

assert sys.version_info >= (3, 12)

def equal(a, b):
    assert s.simplify(a - b) == 0, (a, b)

# Fidelity: Example 3.3, coefficient algebra C^2, frame (p1,p2,p2).
p1, p2 = s.diag(1, 0), s.diag(0, 1)
frame = (p1, p2, p2)
S = sum((x * x.conjugate().T for x in frame), s.zeros(2))
assert S == s.diag(1, 2)
coefficients = tuple(x * (S.inv() * x).conjugate().T for x in frame)
assert coefficients == (p1, p2 / 2, p2 / 2)
assert sum((p - h for p, h in zip(frame, coefficients)), s.zeros(2)) == p2
assert sum((p - p for p in (p1, p2)), s.zeros(2)) == s.zeros(2)
print('FIDELITY Example 3.3: S=diag(1,2); coefficients=(p1,p2/2,p2/2); excess=p2; Riesz excess=0')

# Exact symbolic identities for the continuous-parameter family (t,1).
t, a, b = s.symbols('t a b', real=True)
q = 1 + t**2
T, analysis = s.Matrix([[t, 1]]), s.Matrix([[t], [1]])
assert T * analysis == s.Matrix([[q]])
k = s.Matrix([1, -t])
assert T * k == s.zeros(1, 1)
assert T * s.Matrix([a, b]) == s.Matrix([t*a+b])
# Substituting the equation b=-ta gives every kernel vector as a*(1,-t).
assert s.Matrix([a, b]).subs(b, -t*a) == a*k
P = s.eye(2) - analysis * T / q
assert s.simplify(P - k*k.T/q) == s.zeros(2)
assert s.simplify(P*P - P) == s.zeros(2)
assert s.simplify(T*P) == s.zeros(1, 2)
assert P.T == P
assert s.simplify(s.trace(P)) == 1
h1, h2 = t**2/q, 1/q
equal(h1+h2, 1)
carrier = s.symbols('carrier')
equal((carrier-h1)+(1-h2), carrier)
# Frame bounds and all complex x=p+ir are checked symbolically.
p, r = s.symbols('p r', real=True)
energy = q*(p*p+r*r)
equal(energy-(p*p+r*r), t*t*(p*p+r*r))
equal(2*(p*p+r*r)-energy, (1-t*t)*(p*p+r*r))
print('FAMILY identities: frame multiplier=1+t^2; kernel=A*(1,-t); kernel projection idempotent; diagonal sum=1; excess=carrier(h1)')

# Carrier spectral approximants: n*h/(1+n*h), with exact error and monotonicity.
n = s.symbols('n', integer=True, positive=True)
u = s.symbols('u', real=True, nonnegative=True)
h = u/(1+u)
rn = n*h/(1+n*h)
equal(rn, n*u/(1+(n+1)*u))
equal(1-rn, (1+u)/(1+(n+1)*u))
equal(rn.subs(n, n+1)-rn, u*(1+u)/((1+(n+1)*u)*(1+(n+2)*u)))
assert rn.subs(u, 0) == 0
# Formal symbolic limit for every strictly positive u, separate from u=0.
v = s.symbols('v', positive=True)
assert s.limit(rn.subs(u, v), n, s.oo) == 1
print('CARRIER identities: r_n=n*u/(1+(n+1)*u); r_(n+1)-r_n>=0; error=(1+u)/(1+(n+1)*u); limit=0 at u=0 and 1 at u>0')

# Rational sample checks, including the exact discontinuity sequence t_m=1/m.
points = sorted({Q(j, 32) for j in range(-32, 33)})
for x in points:
    den = 1+x*x
    hs = (x*x/den, 1/den)
    assert Q(1) <= den <= Q(2)
    assert sum(hs) == 1 and hs[1] > 0
    excess = int(x != 0) - hs[0] + 1 - hs[1]
    assert excess == int(x != 0) and 0 <= excess <= 1
    previous = Q(0)
    for N in range(1, 65):
        value = N*hs[0]/(1+N*hs[0])
        assert previous <= value <= int(x != 0)
        previous = value
for m in range(1, 129):
    x = Q(1, m)
    assert x*x/(1+x*x) == Q(1, m*m+1)
    assert 1/(1+x*x) == Q(m*m, m*m+1)
    assert int(x != 0) == 1
print('RATIONAL checks: 65 interval points x 64 approximants; 128 exact t_m=1/m values; sharp bounds attained at t=0,1')

# Arbitrary complex scalar f=u+iv: the same algebra depends only on |f|^2.
z, w = s.symbols('z w', real=True)
f, fc = z+s.I*w, z-s.I*w
modsq = s.expand(f*fc)
assert modsq == z*z+w*w
Tf, af = s.Matrix([[f, 1]]), s.Matrix([[fc], [1]])
kf = s.Matrix([1, -f])
assert s.simplify(Tf*af-s.Matrix([[1+modsq]])) == s.zeros(1)
assert s.simplify(Tf*kf) == s.zeros(1)
equal(modsq/(1+modsq)+1/(1+modsq), 1)
print('GENERAL complex-family identities: (f,1), multiplier=1+|f|^2, kernel generator=(1,-f), coefficient sum=1')
print('ANALYTIC boundary: exact algebra does not certify carrier existence in the bidual, normal extension, continuity contradiction, or all real parameters')
print('ALL_CHECKS_OK')
