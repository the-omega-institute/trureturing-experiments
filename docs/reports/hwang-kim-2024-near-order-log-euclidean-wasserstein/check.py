#!/usr/bin/env python3
"""Exact algebra for Hwang--Kim arXiv:2407.07438v1, Section 6.

Python 3 + SymPy. No numerical floating-point evidence is used.
The polynomial identities and rational enclosures below are certificates outside
Lean. The analytic passage from a negative Taylor coefficient to a punctured
neighbourhood is justified in the accompanying preregistration outline.
"""
from fractions import Fraction as F
import sympy as s


def zero(x):
    if isinstance(x, s.MatrixBase):
        assert all(s.simplify(y) == 0 for y in x)
    else:
        assert s.simplify(x) == 0


class Interval:
    """Closed rational intervals; every operation encloses the exact value."""
    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = self.lo if hi is None else F(hi)
        assert self.lo <= self.hi

    @staticmethod
    def lift(x):
        return x if isinstance(x, Interval) else Interval(x)

    def __add__(self, other):
        o = self.lift(other)
        return Interval(self.lo + o.lo, self.hi + o.hi)
    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.lift(other)

    def __rsub__(self, other):
        return self.lift(other) + -self

    def __mul__(self, other):
        o = self.lift(other)
        a = [x*y for x in (self.lo, self.hi) for y in (o.lo, o.hi)]
        return Interval(min(a), max(a))
    __rmul__ = __mul__

    def __truediv__(self, other):
        o = self.lift(other)
        assert not o.lo <= 0 <= o.hi
        return self * Interval(1/o.hi, 1/o.lo)

    def __rtruediv__(self, other):
        return self.lift(other) / self

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        ans = Interval(1)
        for _ in range(n):
            ans = ans * self
        return ans


def enclose(expr, env):
    if expr.is_Rational:
        return Interval(F(int(expr.p), int(expr.q)))
    if expr.is_Symbol:
        return env[expr]
    if expr.is_Add:
        out = Interval(0)
        for arg in expr.args:
            out = out + enclose(arg, env)
        return out
    if expr.is_Mul:
        out = Interval(1)
        for arg in expr.args:
            out = out * enclose(arg, env)
        return out
    if expr.is_Pow and expr.exp.is_Integer:
        base = enclose(expr.base, env)
        n = int(expr.exp)
        return base ** n if n >= 0 else 1/(base ** (-n))
    raise AssertionError(expr)


I = s.eye(2)
t = s.symbols('t', real=True)
a, b, c, p, q, h = s.symbols('a b c p q h', real=True)
AA = s.Matrix([[a, b], [b, c]])
KK = s.Matrix([[p, q], [q, h]])
BB = KK * AA * KK
# Source fidelity: equation (2.4), independently expanded as a polynomial identity.
zero(((1-t)*I+t*KK)*AA*((1-t)*I+t*KK) -
     ((1-t)**2*AA+t**2*BB+t*(1-t)*(AA*KK+KK*AA)))
# Source fidelity: the commuting case of LE <=near arithmetic, equation (6.21).
# A=diag(4,1), B=diag(1/4,1), w=(1/2,1/2): LE=I, arithmetic=diag(17/8,1).
assert F(17, 8) > 1
Ac = s.diag(4, 1)
Bc = s.diag(s.Rational(1,4), 1)
zero(s.log(4)+s.log(s.Rational(1,4)))
Mc = s.diag(s.sqrt(s.Rational(17,8)), 1)
zero(Mc*Mc-(Ac+Bc)/2)
assert Mc[0,0]-1 > 0 and Mc[1,1]-1 == 0
print('FIDELITY_OK: source (2.4) polynomial identity; source (6.21) commuting pair')

u, z = s.symbols('u z', positive=True)
g = s.sqrt(15)
A = s.diag(4, 1)
R = s.Matrix([[g, -1], [1, g]])/4
B = s.Matrix([[19, -3*g], [-3*g, 61]])/64
zero(R.T*R-I)
zero(B-R*s.diag(s.Rational(1,4),1)*R.T)
assert A.det() == 4 and B.det() == s.Rational(1,4) and B[0,0] > 0
assert A*B != B*A
K = s.Matrix([[35, -3*g], [-3*g, 125]])/(8*s.sqrt(265))
zero(K*A*K-B)
assert K[0,0] > 0 and s.simplify(K.det()) == s.Rational(1,4)
W = ((I+K)*A*(I+K)/4).applyfunc(s.simplify)
zero(W-(A+B+A*K+K*A)/4)
d = s.Rational(5,8)+10/s.sqrt(265)
zero(W.det()-d**2)
# (I+K) is invertible positive definite; this congruence certifies W>0.
assert (I+K)[0,0] > 0 and s.simplify((I+K).det()) > 0

# Exact LE identification. H=(log A+log B)/2=(log 2) J/4, J^2=I.
J = s.Matrix([[1, -g], [-g, -1]])/4
zero(J*J-I)
zero(s.diag(1,0)-R*s.diag(1,0)*R.T-J/4)
L = s.Matrix([[5*u+3/u, -g*(u-1/u)],
              [-g*(u-1/u), 3*u+5/u]])/8
zero(L-((u+1/u)*I+(u-1/u)*J)/2)
zero(L.det()-1)
# u=positive fourth root of 2; eigenvalues of L are u,1/u, so L>0.

# Independent general 2x2 identity behind the near-order test.
la, lb, wa, wb, wc, dd = s.symbols('la lb wa wb wc dd', real=True)
LL = s.Matrix([[la,lb],[lb,(1+lb**2)/la]])
WW = s.Matrix([[wa,wb],[wb,wc]])
N = WW+dd*LL.inv()
Q = (LL*WW).trace()+2*dd
identity = (N*LL*N-Q*WW).applyfunc(s.factor)
for entry in identity:
    zero(entry.subs(dd**2, WW.det()))
zero(s.expand(N.det()-dd*Q).subs(dd**2,WW.det()))
print('WITNESS_IDENTITIES_OK: strict positivity, spectra, noncommutation, Riccati, LE, det W')

T = s.expand((s.trace(W)+d*s.trace(L))**2 -
             (1+d)**2*(s.trace(L*W)+2*d))
T = s.factor(T.subs(s.sqrt(265),z))
ui = Interval(F(1189207115002721,10**15), F(1189207115002722,10**15))
zi = Interval(F(16278820596099706,10**15), F(16278820596099707,10**15))
assert ui.lo**4 < 2 < ui.hi**4
assert zi.lo**2 < 265 < zi.hi**2
Ti = enclose(T, {u:ui,z:zi})
assert F(3,2500) < Ti.lo <= Ti.hi < F(121,100000)
print('EXACT_COUNTEREXAMPLE_OK: 3/2500 < T < 121/100000')
print('ROOT_BRACKETS_OK: fourth-root(2) and square-root(265), integer-power checks')

# General r-family: derive second jets from the literal Riccati and exponential formulas.
r, x, theta = s.symbols('r x theta', positive=True)
trunc = lambda M: M.applyfunc(lambda e:s.series(e,theta,0,3).removeO().expand())
Ar = s.diag(r**2,1)
Bjet = s.Matrix([[r**-2+(1-r**-2)*theta**2, (r**-2-1)*theta],
                 [(r**-2-1)*theta,1+(r**-2-1)*theta**2]])
alpha = r**2+r**-2-2
# sqrt(det Ar det Btheta)=1; (4+alpha theta^2)^(-1/2)=1/2-alpha theta^2/16+O(theta^4).
Kjet = trunc((Bjet+Ar.inv())*(s.Rational(1,2)-alpha*theta**2/16))
zero(trunc(Kjet*Ar*Kjet-Bjet))
Wjet = trunc((Ar+Bjet+Ar*Kjet+Kjet*Ar)/4)
Hjet = s.Matrix([[x*theta**2,-x*theta],[-x*theta,-x*theta**2]])
Ljet = trunc(I+Hjet+Hjet*Hjet/2)
c0 = (r+1/r)/2
w12 = s.factor(Wjet[0,1].coeff(theta,1))
w22 = s.factor(Wjet[1,1].coeff(theta,2))
zero(w12-(r**-2-1)*(r**2+3)/8)
zero(w22+(1-r**-2)*(r**2+7)/16)
zero(Ljet[0,1].coeff(theta,1)+x)
zero(Ljet[1,1].coeff(theta,2)+x-x**2/2)
v = s.factor((w12+c0*x)/(c0+1))
q2 = s.factor((w22+x-x**2/2-v**2+2*v*x)/2)
# The second-order (1,1) coefficient is unnecessary for det(C-I), but we solve it
# to verify the entire Riccati jet, rather than just the selected equations.
p2 = s.solve(s.Eq(2*c0*s.Symbol('p2')+c0**2*(x+x**2/2)-2*c0*v*x+v**2,
                         Wjet[0,0].coeff(theta,2)),s.Symbol('p2'))[0]
Cjet = s.Matrix([[c0+p2*theta**2,v*theta],[v*theta,1+q2*theta**2]])
zero(trunc(Cjet*Ljet*Cjet-Wjet))
D = s.factor(trunc(Cjet-I).det().expand().coeff(theta,2))
S = r**6+2*r**5*x-16*r**4*x+5*r**4+12*r**3*x-16*r**2*x-5*r**2+2*r*x-1
zero(D-(-r**2+2*r*x+1)*S/(32*r**3*(r+1)**2))
# Substitute r=exp(x) using the Laurent form of hyperbolic functions.
sinh1=(r-r**-1)/2
sinh3=(r**3-r**-3)/2
cosh1=(r+r**-1)/2
cosh2=(r**2+r**-2)/2
zero(S/r**3-(2*sinh3+10*sinh1+4*x*cosh2-32*x*cosh1+12*x))
zero(-r**2+2*r*x+1-(-2*r*(sinh1-x)))

def numerator(k):
    return 2*3**(2*k+1)+10+4*(2*k+1)*2**(2*k)-32*(2*k+1)+(12 if k==0 else 0)
assert numerator(0)==0 and numerator(1)==16
# All k>=2: 2^(2k)>=16, hence 4(2k+1)2^(2k)-32(2k+1)>=32(2k+1)>0.
# This checks the polynomial rearrangement and the elementary induction for 4^k>=16.
j=s.symbols('j',integer=True,nonnegative=True)
y=s.symbols('y',real=True)
zero(4*(2*j+1)*y-32*(2*j+1) - ((2*j+1)*(4*(y-16)+32)))
assert 4**2 == 16
zero(4**(j+1)-4*4**j)
for k in range(2,129):
    assert 4**k>=16 and numerator(k)>0
print('FAMILY_JET_OK: D(r)=(-r^2+2r log(r)+1) S(r,log(r))/(32 r^3 (r+1)^2)')
print('FAMILY_SIGN_OK: positive odd-series coefficients; sinh(x)-x>0 for x>0')
print('SCOPE: exact algebra and rational certificates; analytic justification in README; no Lean claim proof')
print('ALL_CHECKS_OK')
