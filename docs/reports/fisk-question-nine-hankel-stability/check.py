#!/usr/bin/env python3
"""Exact checks for Fisk, arXiv:0808.1850v1, Question 9 and matrix (4)."""
from itertools import permutations
from fractions import Fraction
import sympy as S

x, y, A, B, w, r, t = S.symbols('x y A B w r t')
a = S.symbols('a', positive=True)
b, c = S.symbols('b c', positive=True)


def zero(expr):
    assert S.cancel(S.expand(expr)) == 0, expr


def determinant(rows):
    """Leibniz expansion, independent of SymPy's determinant implementation."""
    n = len(rows)
    ans = 0
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        ans += (-1) ** inversions * S.prod(rows[i][p[i]] for i in range(n))
    return S.expand(ans)


# Fidelity: Remark 1's degree-two answer to Question 1 and its discriminant.
coeffs = [b*c, b+c, S.Integer(1)]
def ci(i):
    return coeffs[i] if 0 <= i < len(coeffs) else S.Integer(0)
F = sum(S.det(S.Matrix([[ci(i), ci(i+1)], [ci(i-1), ci(i)]]))*x**i
        for i in range(3))
zero(F - (b*b*c*c + (b*b+b*c+c*c)*x + x*x))
zero(S.discriminant(F, x) - (b*b-b*c+c*c)*(b*b+3*b*c+c*c))
print('FIDELITY_OK: source Remark 1, degree-two transform and discriminant (symbolic b,c>0)')

# A concrete instance of Remark 8's d=1 assertion.
f_repeated = S.Poly((y+x+1)**4, y)
h2 = determinant([[f_repeated.nth(0), f_repeated.nth(1)],
                  [f_repeated.nth(1), f_repeated.nth(2)]])
zero(h2 + 10*(x+1)**6)
print('FIDELITY_OK: Remark 8, f=(x+y+1)^4, i=0,d=1: F=-10(x+1)^6')

# The full triangular input support is positive for every real a>0.
f = S.Poly(S.expand((a*x+y+1)**2*(x+y+a)**2), x, y)
assert set(f.monoms()) == {(i,j) for i in range(5) for j in range(5-i)}
for _, coeff in f.terms():
    assert all(v > 0 for v in S.Poly(coeff, a).all_coeffs() if v != 0)
xr, yr = S.symbols('xr yr', real=True)
xi, yi = S.symbols('xi yi', positive=True)
assert S.im(a*(xr+S.I*xi)+(yr+S.I*yi)+1) == a*xi+yi
assert (a*xi+yi).is_positive
assert S.im((xr+S.I*xi)+(yr+S.I*yi)+a) == xi+yi
assert (xi+yi).is_positive
print('INPUT_OK: 15 positive support coefficients and positive imaginary parts of both affine factors')

fy = S.Poly((y+A)**2*(y+B)**2, y)
cs = [fy.nth(j) for j in range(7)]
rows = [[cs[i+j] for j in range(4)] for i in range(4)]
H = determinant(rows)
zero(H - S.det(S.Matrix(rows)))
h = 5*A**4+8*A**3*B+9*A**2*B**2+8*A*B**3+5*B**4
zero(H-h)
zero(H-(cs[3]**4-3*cs[2]*cs[3]**2+cs[2]**2+2*cs[1]*cs[3]-cs[0]))
Ha = S.expand(h.subs({A:a*x+1, B:x+a}))
assert S.Poly(Ha.subs(x,0), a).all_coeffs() == [5,8,9,8,5]
print('DETERMINANT_OK: prescribed i=0,d=3 four-by-four determinant equals 5A^4+8A^3B+9A^2B^2+8AB^3+5B^4')

q = 5*w**4+8*w**3+9*w*w+8*w+5
zero(q-w*w*(5*(w+1/w)**2+8*(w+1/w)-1))
rp = (S.sqrt(21)-4)/10
rm = (-S.sqrt(21)-4)/10
assert bool(-1 < rm) and bool(rm < 0) and bool(0 < rp) and bool(rp < 1)
zero(q-5*(w*w-2*rp*w+1)*(w*w-2*rm*w+1))
for rr in (rp,rm):
    zero((rr+S.I*S.sqrt(1-rr*rr))**2-2*rr*(rr+S.I*S.sqrt(1-rr*rr))+1)
print('AUXILIARY_ROOTS_OK: q factors into two distinct unit-circle conjugate pairs')

# Independent real-quadratic factorization certifies the full parameter classification.
def K(rr):
    return S.expand((a*a+1-2*rr*a)*(x*x+1)+(4*a-2*rr*(a*a+1))*x)
zero(Ha-5*K(rp)*K(rm))
zero(S.discriminant(K(r),x)-4*(r*r-1)*(a*a-1)**2)
zero(a*a+1-2*r*a-((a-r)**2+1-r*r))
zero(S.resultant(K(rp),K(rm),x)-4*(rp-rm)**2*(a*a-1)**4)
# For a != 1, each quadratic has a nonreal pair with real part
# (r*(a*a+1)-2*a)/(a*a+1-2*r*a), and its denominator is positive.
threshold = 16+4*S.sqrt(21)
zero(2/rp-threshold)
zero(rp*(a*a+1)-2*a-a*rp*(a+1/a-threshold))
zero(Ha.subs(a,1)-35*(x+1)**4)
assert bool(threshold < 35)
s = S.symbols('s', nonnegative=True)
positive_numerator = S.Poly(S.expand((rp*(a*a+1)-2*a).subs(a,s+35)),s)
assert all(bool(v > 0) for v in positive_numerator.all_coeffs())
assert all(bool(v < 0) for v in S.Poly(rm*(a*a+1)-2*a,a).all_coeffs())
zero(Ha-a**4*Ha.subs(a,1/a))
print('CLASSIFICATION_OK: a>0, stable iff a+1/a <= 16+4*sqrt(21); a=1 treated separately')
print('FAMILY_OK: all a>=35 and all 0<a<=1/35 have exactly two simple right-half-plane roots')

# Verify the fractional-linear transport without evaluating approximate roots.
z = (a*w-1)/(a-w)
zero((a-w)**4*Ha.subs(x,z)-(a*a-1)**4*q)
u,v = S.symbols('u v', real=True)
D = (a-u)**2+v*v
N = S.expand((a*(u+S.I*v)-1)*(a-u+S.I*v))
zero(S.re(N)-((a*a+1)*u-2*a)-a*(1-u*u-v*v))
zero(S.im(N)-(a*a-1)*v)
zero(S.expand((a*u-1)**2+(a*v)**2-D)-(a*a-1)*(u*u+v*v-1))
print('TRANSPORT_OK: determinant, real-part and unit-circle identities are exact polynomial identities')

# Explicit a=35 algebraic root, reduced by t^2=1-r^2 then r=rp.
L = a*a+1-2*r*a
root_numerator = r*(a*a+1)-2*a+S.I*(a*a-1)*t
residual = S.expand(L**2*K(r).subs(x,root_numerator/L))
remainder = S.rem(residual,t*t-(1-r*r),t)
zero(remainder)
L35 = L.subs({a:35,r:rp})
R35 = (r*(a*a+1)-2*a).subs({a:35,r:rp})
assert bool(L35>0) and bool(R35>0)
print('EXACT_WITNESS_OK: a=35, Re(z) =',S.radsimp(R35/L35),'> 0')

for n in range(1,51):
    fs = S.Poly(f.as_expr().subs(a,n), y)
    direct = determinant([[fs.nth(i+j) for j in range(4)] for i in range(4)])
    zero(direct-Ha.subs(a,n))
for value in (Fraction(1,35),Fraction(1,2),Fraction(1),Fraction(2),Fraction(35)):
    aa=S.Rational(value.numerator,value.denominator)
    assert bool((rp*(aa*aa+1)-2*aa)>0) == (aa in [S.Rational(1,35),S.Integer(35)])
print('SMALL_CASES_OK: 50 integer determinant recomputations and five exact rational root-sign checks')
print('ALL_CHECKS_OK')
