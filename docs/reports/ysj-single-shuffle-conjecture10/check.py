#!/usr/bin/env python3
"""Exact checks for Yun--Sra--Jadbabaie, arXiv:2103.07079v1, (25)--(26).

The three projections and near-identity construction are from Appendix A of
that source and Peng, arXiv:2607.22620v1, Section 2. The polynomial gap is Peng's
violation lemma (Lemma 6). The fixed-determinant expectation-of-norm reduction is the candidate
being checked. This is exact symbolic/rational evidence, not a Lean proof.
Requires Python 3 and SymPy. No random inputs or floating-point arithmetic.
"""
from fractions import Fraction
from itertools import permutations
from math import isqrt
import sympy as S

q = S.Symbol('q', real=True)
I = S.eye(2)
u = [S.Matrix([1, 0]), S.Matrix([S.Rational(1, 2), S.sqrt(3)/2]),
     S.Matrix([-S.Rational(1, 2), S.sqrt(3)/2])]
P = [v*v.T for v in u]
perms = list(permutations(range(3)))


def eq(x, y):
    if isinstance(x, S.MatrixBase):
        assert x.shape == y.shape
        assert all(S.expand(z) == 0 for z in x-y), (x, y)
    else:
        assert S.expand(x-y) == 0, S.factor(x-y)


def word(factors, order):
    result = I
    for i in order:
        result = result*factors[i]
    return result.applyfunc(S.expand)


def frob2(x):
    return S.expand(S.trace(x.T*x))


for p in P:
    eq(p*p, p)
    eq(p.T, p)
    eq(p.trace(), 1)
print('PROJECTIONS: symmetry, idempotence, trace 1 checked exactly')

# Fidelity: reproduce the displayed Appendix A mean-of-powers identity,
# not either target conjecture. The source labels its projectors P2,P3,P1.
source_matrices = [P[1], P[2], P[0]]
source_words = [word(source_matrices, p) for p in perms]
for k in range(1, 9):
    eq(sum((t**k for t in source_words), S.zeros(2))/6,
       S.Rational(-1, 8)**k * I/2)
print('SOURCE FIDELITY: Appendix A, mean of powers = (-1/8)^K I/2, K=1,...,8')

B = [q*I+(1-q)*p for p in P]
for m in B:
    eq(m.T, m)
    eq(m.det(), q)
    eq(m.trace(), 1+q)
    eq(m*m-(1+q)*m+q*I, S.zeros(2))
T = [word(B, p) for p in perms]
d0 = q**3
e = (1-q)**6/32
c = d0+e
t = (q+1)*(-q**2+10*q-1)/8
a = d0**2+t**2*e
b = c**2
D = q**6
for x in T:
    eq(x.det(), d0)
    eq(x.trace(), t)
    eq(frob2(x), 2*c)
    eq(x*x-t*x+d0*I, S.zeros(2))
    eq(frob2(x*x), 2*a)
    eq((x*x).det(), D)
print('ONE EPOCH: all 6 determinant, trace, Frobenius, Cayley-Hamilton identities')
print('REPEATED EPOCHS: all 6 squared-product Frobenius norms = 2a; determinant q^6')

eq(sum((x.T*x for x in T), S.zeros(2))/6, c*I)
eq(sum((x*x.T for x in T), S.zeros(2))/6, c*I)
products = [x*y for x in T for y in T]
for x in products:
    eq(x.det(), D)
eq(sum(frob2(x) for x in products)/36, 2*b)
print('INDEPENDENT EPOCHS: all 36 determinants q^6; mean squared Frobenius norm = 2b')

g = 42*q**2-q**4-4*q**3-4*q-1
gap = (1-q)**6*(q+1)**2*g/2048
eq(a-b, gap)
# The four displayed summands are nonnegative for 1/2 <= q <= 1.
eq(g-25*q**2, q**2*(1-q**2)+4*q**2*(1-q)
   +4*q*(2*q-1)+(2*q-1)*(2*q+1))
assert S.Rational(1, 2)**2 * 25 > 0
assert S.factor(g.subs(q, 1)) == 32
assert (a-b).subs(q, S.Rational(3, 4)) == S.Rational(221039, 34359738368)
print('FAMILY GAP: a-b = (1-q)^6 (q+1)^2 g(q)/2048')
print('SIGN CERTIFICATE: g(q) >= 25q^2 > 0 on [1/2,1]; gap > 0 for q < 1')
print('RATIONAL SAMPLE: q=3/4; a-b=221039/34359738368')

# Symbolic algebra underlying fixed-determinant singular-value transfer.
al, be = S.symbols('alpha beta', positive=True)
eq((al+be)**2, al**2+be**2+2*al*be)
eq((al-be)**2, al**2+be**2-2*al*be)
s, dd = S.symbols('s D', positive=True)
H = (S.sqrt(s+2*dd)+S.sqrt(s-2*dd))/2
F = (s+S.sqrt((s-2*dd)*(s+2*dd)))/2
# These are identities on s > 2D, where the displayed square roots are positive.
eq(S.diff(H, s, 2), -1/(8*(s+2*dd)**S.Rational(3, 2))
   -1/(8*(s-2*dd)**S.Rational(3, 2)))
r = S.Symbol('r', positive=True)
F_interior = (2*dd+r+S.sqrt(r*(r+4*dd)))/2
assert S.simplify(S.diff(F_interior, r, 2)
                  + 2*dd**2/(r*(r+4*dd))**S.Rational(3, 2)) == 0
print('TRANSFER: exact singular-value square identities and concavity derivatives checked')

# Direct checks of spectral-norm expectations at q=3/4, with rational
# enclosures of every square root. No numerical sign oracle is used.
SCALE = 10**40


def rational(x):
    x = S.cancel(x)
    assert x.is_Rational
    return Fraction(int(x.p), int(x.q))


def sqrt_bounds(x):
    assert x >= 0
    k = isqrt((x.numerator*SCALE*SCALE)//x.denominator)
    lo, hi = Fraction(k, SCALE), Fraction(k+1, SCALE)
    assert lo*lo <= x <= hi*hi
    return lo, hi


def norm_bounds(sf, det):
    l1, h1 = sqrt_bounds(sf+2*det)
    l2, h2 = sqrt_bounds(sf-2*det)
    return (l1+l2)/2, (h1+h2)/2


def norm_square_bounds(sf, det):
    lo, hi = sqrt_bounds((sf-2*det)*(sf+2*det))
    return (sf+lo)/2, (sf+hi)/2


sample = S.Rational(3, 4)
det = rational(D.subs(q, sample))
srep = rational((2*a).subs(q, sample))
sind = [rational(frob2(x).subs(q, sample)) for x in products]
assert all(z >= 2*det for z in sind+[srep])
for label, bound in [('(25)', norm_bounds), ('(26)', norm_square_bounds)]:
    repeated_lo, _ = bound(srep, det)
    independent_hi = sum(bound(z, det)[1] for z in sind)/36
    lower_gap = repeated_lo-independent_hi
    assert lower_gap > Fraction(1, 10**8), (label, lower_gap)
    print('DIRECT SPECTRAL CHECK '+label+': repeated expectation - independent expectation > 1/10^8')

# Tolerance selection is by two exhaustive exact branches, not a finite grid.
eta = S.Symbol('eta', real=True)
q_small = 1-eta/2
q_large = S.Rational(3, 4)
eq(q_small-(1-eta), eta/2)
eq(q_large-(1-eta), eta-S.Rational(1, 4))
print('TOLERANCE: q=1-min(eta/2,1/4); branches eta<=1/2 and eta>=1/2')
print('For every 0<eta<=1, 3/4<=q<1 and q>=1-eta; B_i have eigenvalues q,1')
print('ALL_CHECKS_OK')
