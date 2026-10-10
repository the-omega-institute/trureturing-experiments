# Exact check of the scalar inequality that arXiv:2206.06653v1 Section 1 quotes as the Kushel–Tyaglov theorem (and that Conjecture 2.4 lifts),
# at d = 3 with roots a = (1+i, 1-i, 1) in the complex numbers.
import sympy as sp
I=sp.I; z=sp.symbols('z')
a=[1+I,1-I,sp.Integer(1)]; d=3
P=sp.expand(sp.prod([z-x for x in a]))
dP=sp.diff(P,z)
b=sp.solve(sp.Eq(dP,0),z)
assert sp.expand(dP-3*(z-b[0])*(z-b[1]))==0
ab2=lambda x: sp.simplify(x*sp.conjugate(x))
S=sum(a); T=sum(x**2 for x in a)-S**2/d**2
lhs=sp.nsimplify(sum(ab2(x)**2 for x in b))
rhs=sp.nsimplify(sp.Rational(d-6,d)*sum(ab2(x)**2 for x in a)+sp.Rational(1,d**2)*sum(ab2(x) for x in a)**2+sp.Rational(1,d**2)*ab2(T)
     +sp.Rational(2,d)*sum(ab2(x)*ab2(x+S/d) for x in a)-sp.Rational(4,d**3)*sum(ab2(x)*ab2(S) for x in a))
print('P =',P,'; zeros of P\' =',b)
print('lhs =',lhs,' rhs =',rhs,' rhs - lhs =',sp.simplify(rhs-lhs))
assert sp.simplify(rhs-lhs)==sp.Rational(-4,9)
print('OK')
