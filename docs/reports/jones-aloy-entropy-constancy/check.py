#!/usr/bin/env python3
"""Exact checks for the entropy-constancy candidate in arXiv:2506.14526.

Rational logarithm enclosures use the convergent atanh series and an explicit
geometric tail bound. SymPy checks the identities used for the entire parameter
family. These checks support the mathematical argument in README.md; they are
not a formal proof of an optimization theorem or a numerical optimizer.
"""
from fractions import Fraction as Q
import ast
from pathlib import Path
import sympy as S


def same(a, b):
    assert S.simplify(a - b) == 0, (a, b)


def log_reduced(x, terms=40):
    """Enclose ln(x) for 1 <= x <= 2, using rational arithmetic only."""
    assert Q(1) <= x <= Q(2)
    z = (x - 1) / (x + 1)
    lo = 2 * sum((z ** (2*k + 1) / (2*k + 1)
                  for k in range(terms)), Q(0))
    tail = 2 * z ** (2*terms + 1) / ((2*terms + 1) * (1 - z*z))
    return lo, lo + tail


def log_interval(x):
    assert x > 0
    k = 0
    while x < 1:
        x *= 2
        k -= 1
    while x > 2:
        x /= 2
        k += 1
    lo, hi = log_reduced(x)
    l2, u2 = log_reduced(Q(2))
    return (lo + k*l2, hi + k*u2) if k >= 0 else (lo + k*u2, hi + k*l2)


def entropy_nats_interval(p):
    assert 0 <= p <= 1
    if p in (0, 1):
        return Q(0), Q(0)
    lp, up = log_interval(p)
    lr, ur = log_interval(1-p)
    return -p*up-(1-p)*ur, -p*lp-(1-p)*lr


def coarse_lower(x, denominator=10**9):
    return Q((x.numerator * denominator) // x.denominator, denominator)


# Reject accidental floating-point constants in the evidence program itself.
assert not any(isinstance(x, ast.Constant) and isinstance(x.value, float)
               for x in ast.walk(ast.parse(Path(__file__).read_text())))

# Fidelity: source Appendix E(ii), Lemma 3 midpoint identity and the reflection
# invariance from Lemma 2, reproduced at every rational correlation, symbolically.
x, y = S.symbols('x y', real=True)
C = S.Matrix([x, y])
R = S.Matrix([-y, -x])
M = S.Matrix([(x-y)/2, -(x-y)/2])
assert S.simplify((C+R)/2-M) == S.zeros(2, 1)
same(R[0]-R[1], x-y)
h = lambda z: -((1+z)/2)*S.log((1+z)/2) - ((1-z)/2)*S.log((1-z)/2)
same(h(-x), h(x))
fid = lambda u,v: (S.sqrt(1+u)*S.sqrt(1+v)+S.sqrt(1-u)*S.sqrt(1-v))/2
same(fid(-y,-x), fid(x,y))
print('FIDELITY_OK: source Lemmas 2-3 reflection and midpoint identities (all x,y)')

s, a, L, u, t = S.symbols('s a L u t', real=True)
q = s*s
B = S.Matrix([1, 1-2*q])
mid = S.Matrix([q, -q])
D = S.Matrix([s, -s])
same(B[0]-B[1], mid[0]-mid[1])
ensemble_mean = s*D + (1-s)/2*S.Matrix([1,1]) + (1-s)/2*S.Matrix([-1,-1])
assert S.simplify(ensemble_mean-mid) == S.zeros(2,1)
E = S.symbols('E', real=True)
same(s*E+(1-s)/2*E+(1-s)/2*E, E)
same(1-q*q-(1-q), q*(1-q))
print('FAMILY_ENSEMBLE_OK: normalized nonnegative weights on 24/25 <= s < 1; mean M and exact energy E')

# Supporting tangent: f(p)=asin(sqrt(p)). On p>1/2 its second derivative is
# positive; on p<1/2 it is negative. The endpoint certificate is uniform in s.
p = S.symbols('p', positive=True)
f = S.asin(S.sqrt(p))
same(S.diff(f,p), 1/(2*S.sqrt(p)*S.sqrt(1-p)))
same(S.diff(f,p,2), (2*p-1)/(4*p**S.Rational(3,2)*(1-p)**S.Rational(3,2)))
same(49*s*s-576*(1-s*s), 625*s*s-576)
endpoint_poly = S.Poly(S.expand((625*s*s-576).subs(s,u+S.Rational(24,25))),u)
assert endpoint_poly.all_coeffs() == [625,1200,0]
assert Q(22,7) < Q(24,7)  # with the classical strict bound pi < 22/7
print('TANGENT_OK: exact derivatives and 625(s-24/25)^2+1200(s-24/25) >= 0')

# Entropy-gap algebra. L=log2(1/a), 0<a<=1/50. The integer inequality proves
# log2(50)>11/2. The atanh series proves ln(2)>2/3 exactly.
assert 50**2 > 2**11
log2_lo, _ = log_interval(Q(2))
assert log2_lo > Q(2,3)
same(1-(1-2*a)**2, 4*a*(1-a))
lower = 2*a*(1-a)*(L-2)
upper = (1-2*a)*a*(L+S.Rational(3,2))
same(lower-upper, a*(L-S.Rational(11,2)+7*a))
same(a*(L-S.Rational(11,2)+7*a)-7*a*a, a*(L-S.Rational(11,2)))
# The scalar log estimate is witnessed by g(0)=0 and g'(a)=-ln(1-a)>0.
g = a+(1-a)*S.log(1-a)
same(S.diff(g,a), -S.log(1-a))
same(g.subs(a,0),0)
print('FAMILY_GAP_OK: gap > a[L-11/2+7a] > 7a^2 > 0 for every 0<a<=1/50')

# Independent rational enclosure of the entropy gap at the exact witness.
s0 = Q(24,25)
q0 = s0*s0
assert (Q(1),1-2*q0) == (Q(1),Q(-527,625))
assert (q0,-q0) == (Q(576,625),Q(-576,625))
assert (1+s0)/2 == Q(49,50)
bq_lo,bq_hi = entropy_nats_interval(q0)
bp_lo,bp_hi = entropy_nats_interval((1+s0)/2)
gap_lo = bq_lo/2-s0*bp_hi
assert gap_lo > 0
assert gap_lo > Q(7,2500)*log_interval(Q(2))[1]
print('EXACT_WITNESS_OK: B=(1,-527/625), M=(576/625,-576/625), c=7/25')
print('EXACT_LOG_GAP_NATS_LOWER_BOUND:',coarse_lower(gap_lo))

# Pauli realization, independent of the entropy computation. The evolved state
# is (c,s), c^2=1-s^2; E sigma_y has variance E^2 in (1,0).
c = S.symbols('c', positive=True)
X=S.Matrix([[0,1],[1,0]])
Y=S.Matrix([[0,-S.I],[S.I,0]])
Z=S.diag(1,-1)
I=S.eye(2)
psi0=S.Matrix([1,0]);psi1=S.Matrix([c,s])
def reduce(expr):
    return S.expand(expr).subs(c**2,1-s*s).expand()
O=-t*s*c*X+(1-t*c*c)*Z
same(reduce((psi1.T*Z*psi1)[0]),1-2*q)
same(reduce((psi1.T*(-c*X+s*Z)*psi1)[0]),-s)
assert S.simplify((E*Y)**2-E*E*I)==S.zeros(2)
same((psi0.T*E*Y*psi0)[0],0)
normsq=1-c*c*t*(2-t)
for i in range(2):
 for j in range(2):same(reduce((O*O-normsq*I)[i,j]),0)
Bt=(1-t)*B+t*mid
same(reduce((psi0.T*O*psi0)[0]),reduce(Bt[0]))
same(reduce((psi1.T*O*psi1)[0]),reduce(Bt[1]))
for sign in (-1,1):same(reduce(((I+sign*O)/2).det()),reduce(c*c*t*(2-t)/4))
print('INTERIOR_REALIZATION_OK: correlations B_t; effect trace=1 and determinant=c^2 t(2-t)/4>0')
print('SCOPE: boundary rigidity and compact-roof lower semicontinuity use the accompanying analytic argument; no numerical optimizer or Lean proof is claimed')
print('ALL_CHECKS_OK')
