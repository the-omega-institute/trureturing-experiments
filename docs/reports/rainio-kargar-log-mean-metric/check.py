#!/usr/bin/env python3
"""Exact source-fidelity and proposed-refutation checks; Python 3 + SymPy.

Source: Rainio--Kargar, arXiv:2308.06576v1, (4.1)--(4.2), Lemma 4.5,
Conjecture 4.6. No floating point or numerical logarithm evaluation is used.
"""
from fractions import Fraction
import sympy as S

Q = S.Rational
c, r, ell = S.symbols('c r ell', positive=True)
s = S.symbols('s', positive=True)
u = S.symbols('u', nonnegative=True)
a, b, t = S.symbols('a b t', positive=True)


def same(actual, expected):
    assert S.cancel(S.expand(actual - expected)) == 0, (actual, expected)
    assert not actual.atoms(S.Float)


def positive_rational(expr, variable):
    """Certify positivity for variable > 0 by positive polynomial coefficients."""
    numerator, denominator = S.fraction(S.cancel(expr))
    for part in (numerator, denominator):
        coefficients = S.Poly(part, variable).all_coeffs()
        assert all(coefficient >= 0 for coefficient in coefficients)
        assert any(coefficient > 0 for coefficient in coefficients)


def norm(v):
    return S.sqrt(sum(x*x for x in v))


def distance(v, w):
    return norm(tuple(x-y for x, y in zip(v, w)))


def logarithmic_mean(x, y):
    if x == y:
        return x
    return S.cancel((x-y)/(S.log(x)-S.log(y)))


def raw_distance(length, dx, dy):
    return S.cancel(length/(c*logarithmic_mean(dx, dy)))


def bounded_distance(length, dx, dy):
    return S.cancel(length/(length+2*c*logarithmic_mean(dx, dy)))


# Analytic certificate for the reused logarithmic-mean bound. For 0<a<b,
# a<=t<=b the chord gap is nonnegative. Integrating 1/t on [a,b] gives
# log(b)-log(a) <= (b-a)*(a+b)/(2*a*b). This independently explains the
# inequality used from the existing D5 kernel sandwich, without approximations.
chord = (a+b-t)/(a*b)
same(chord-1/t, (t-a)*(b-t)/(a*b*t))
same(S.diff(S.log(t), t), 1/t)
same(S.integrate(chord, (t, a, b)), (b-a)*(a+b)/(2*a*b))
same((chord-1/t).subs(t, a+(b-a)*u),
     (b-a)**2*u*(1-u)/(a*b*(a+(b-a)*u)))
log_two_upper = Fraction(3, 4)
assert log_two_upper < 1
print('Log certificate: chord minus 1/t = (t-a)*(b-t)/(a*b*t) >= 0 on [a,b].')
print('Consequences: log(2) <= 3/4; -log(1-r) <= r*(2-r)/(2*(1-r)).')

# Fidelity: the unbounded-function assertion of source Lemma 4.5, for every
# c>0. At x=(1/2,0), y=(-1/2,0), z=(0,0), the boundary distances are 1/2,1/2,1.
x, y, z = (Q(1,2), S.Integer(0)), (-Q(1,2), S.Integer(0)), (S.Integer(0), S.Integer(0))
dx, dy, dz = 1-norm(x), 1-norm(y), 1-norm(z)
assert (dx, dy, dz) == (Q(1,2), Q(1,2), 1)
phi_xy = raw_distance(distance(x, y), dx, dy)
phi_xz = raw_distance(distance(x, z), dx, dz)
same(phi_xy, 2/c)
same(phi_xz, S.log(2)/c)
fidelity_defect_bound = (2-2*Q(log_two_upper.numerator, log_two_upper.denominator))/c
same(fidelity_defect_bound, 1/(2*c))
positive_rational(fidelity_defect_bound, c)
print('FIDELITY Lemma 4.5 (unbounded phi): defect = 2*(1-log(2))/c >= 1/(2*c) > 0.')

# Source (4.2), fixed bounded witness c=2; all lengths and depths are exact.
Dxy = bounded_distance(distance(x, y), dx, dy).subs(c, 2)
Dxz = bounded_distance(distance(x, z), dx, dz).subs(c, 2)
same(Dxy, Q(1,3))
same(2*Dxz, 2*S.log(2)/(S.log(2)+4))
q = S.symbols('q', positive=True)
same(S.diff(2*q/(q+4), q), 8/(q+4)**2)
fixed_upper = 2*Q(3,4)/(Q(3,4)+4)
assert fixed_upper == Q(6,19)
assert Dxy-fixed_upper == Q(1,57)
print('FIXED BALL c=2, r=1/2: Dxy=1/3, Dxz+Dzy <= 6/19, defect >= 1/57.')

# Uniform ball family. Conditions: 0<r<1 and ell=-log(1-r)>0.
# Check the literal source formula first, then the rational comparison.
side_depth = 1-r
source_Dxy = bounded_distance(2*r, side_depth, side_depth)
source_Dxz = bounded_distance(r, side_depth, S.Integer(1))
same(source_Dxy, r/(r+c*(1-r)))
same(source_Dxz, (-S.log(1-r))/(-S.log(1-r)+2*c))
ball_xy = r/(r+c*(1-r))
ball_xz = ell/(ell+2*c)
triangle_numerator = 2*c*r-ell*(r+2*c*(1-r))
same(ball_xy-2*ball_xz,
     triangle_numerator/((r+c*(1-r))*(ell+2*c)))
upper_ell = r*(2-r)/(2*(1-r))
threshold = (2-r)/(2*(1-r))
radius = (c-1)/(2*c-1)
same(c-threshold.subs(r, radius), (c-1)*(2*c-1)/(2*c))
for positive_quantity in (radius, 1-radius, Q(1,2)-radius):
    positive_rational(positive_quantity.subs(c, 1+s), s)
# Monotonicity in ell makes this rational positive quantity a defect lower bound.
same(S.diff(2*ell/(ell+2*c), ell), 4*c/(ell+2*c)**2)
ball_lower = S.factor((ball_xy-2*ball_xz.subs(ell, upper_ell)).subs(r, radius))
positive_rational(ball_lower.subs(c, 1+s), s)
print('BALL FAMILY c>1: r=(c-1)/(2*c-1), 0<r<1/2.')
print('BALL FAMILY defect lower bound:', ball_lower)

# Half-space family in dimension 2, extended by zero horizontal coordinates
# in larger dimensions. Last coordinate is the source boundary distance.
xh, yh, zh = (Q(3,4),Q(1,4)), (-Q(3,4),Q(1,4)), (S.Integer(0),S.Integer(1))
assert distance(xh, yh) == Q(3,2)
assert distance(xh, zh) == 3*S.sqrt(2)/4
assert distance(yh, zh) == 3*S.sqrt(2)/4
assert all(v[-1] > 0 for v in (xh,yh,zh))
B = 2*S.sqrt(2)*S.log(2)
half_xy = bounded_distance(distance(xh, yh), xh[-1], yh[-1])
half_xz = bounded_distance(distance(xh, zh), xh[-1], zh[-1])
same(half_xy, 3/(3+c))
same(half_xz, B/(B+2*c))
bsym = S.symbols('B', positive=True)
same(3/(3+c)-2*bsym/(bsym+2*c),
     (6*c-(2*c+3)*bsym)/((3+c)*(bsym+2*c)))
# From log(2)<=3/4: B<=3*sqrt(2)/2<24/11. The strict radical
# comparison is certified by positive quantities and their integer squares.
assert 2 < Q(16,11)**2
assert Q(3,2)*Q(16,11) == Q(24,11)
same(6*c/(2*c+3)-Q(24,11), 18*(c-4)/(11*(2*c+3)))
# An alternate strictly positive affine certificate also covers c=4 exactly.
Bupper = 3*S.sqrt(2)/2
assert 24**2 > Q(33,2)**2*2
assert 6**2 > 3**2*2
same((6*c-(2*c+3)*Bupper).subs(c, 4+s),
     (24-Q(33,2)*S.sqrt(2))+s*(6-3*S.sqrt(2)))
print('HALF-SPACE FAMILY c>=4: Dxy=3/(3+c), Dxz=Dzy=B/(B+2*c).')
print('HALF-SPACE bound: B<=3*sqrt(2)/2<24/11<=6*c/(2*c+3); strict triangle failure.')

# Exact zero-padding checks of source Euclidean lengths in small dimensions.
for dimension in (2,3,4,7):
    def pad(v):
        return (v[0],)+(S.Integer(0),)*(dimension-2)+(v[-1],)
    for v,w in ((xh,yh),(xh,zh),(yh,zh)):
        assert distance(pad(v),pad(w)) == distance(v,w)
    for v,w in ((x,y),(x,z),(y,z)):
        assert distance(pad(v),pad(w)) == distance(v,w)
print('Zero-padding Euclidean-length checks: dimensions 2,3,4,7, exact.')
print('No claim is made here about metricity in the remaining parameter ranges.')
print('ALL_CHECKS_OK')
