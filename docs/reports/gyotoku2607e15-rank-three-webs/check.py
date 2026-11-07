#!/usr/bin/env python3
"""Exact symbolic checks for Gyotoku, arXiv:2607.20646v1, Example 15.

This is an algebra/calculus check, not a Lean proof or a classification theorem.
Every number used as mathematical evidence is exact. Requires SymPy 1.14.0.
"""
import sympy as s


def zero(expr):
    result = s.factor(s.simplify(expr))
    assert result == 0, result


x, y, u, v, lam = s.symbols('x y u v lambda', real=True)
alpha, beta, delta, t = s.symbols('alpha beta delta t', positive=True)
phi, psi = s.Function('phi'), s.Function('psi')
F = s.Matrix([phi(x+y), psi(x)-phi(x+y)])
zero(F.jacobian([x,y]).det() + s.diff(phi(x+y), x)*s.diff(psi(x), x))
print('FIDELITY equation (18): determinant = -phi_prime(x+y)*psi_prime(x)')

# Remark 14's three Kac--Bernstein relations, exactly as printed.
for rel in [x+y-u, x-y-v, 2*x*x+2*y*y-u*u-v*v]:
    zero(rel.subs({u:x+y, v:x-y}, simultaneous=True))
print('FIDELITY Remark 14: all three printed Gaussian relations vanish')

# Example 15's footnote 4, including normalized input/output densities.
U, V = s.symbols('U V', positive=True)
expF = s.Matrix([s.exp(x+y), s.exp(x+1)-s.exp(x+y)])
expInv = {x:s.log(U+V)-1, y:1+s.log(U/(U+V))}
zero(expF[0].subs(expInv, simultaneous=True)-U)
zero(expF[1].subs(expInv, simultaneous=True)-V)
zero(expF.jacobian([x,y]).det()+s.exp(2*x+y+1))
S, Z = s.symbols('S Z', positive=True)
assert s.integrate(S*s.exp(-S), (S,0,s.oo)) == 1
assert s.integrate(s.Integer(1), (Z,0,1)) == 1
zero((S**2*Z*s.exp(-S)).subs({S:U+V,Z:U/(U+V)}, simultaneous=True)
     /(U*(U+V))-s.exp(-U)*s.exp(-V))
print('FIDELITY footnote 4: inverse, Jacobian and normalized beta--gamma density identity')

# Primary counterexample family: lambda>0, phi(t)=t, psi(t)=(lambda+1)t.
G = s.Matrix([x+y, lam*x-y])
J = G.jacobian([x,y])
assert J.det() == -lam-1
inv = {x:(u+v)/(lam+1),y:(lam*u-v)/(lam+1)}
for a,b in zip(G.subs(inv, simultaneous=True),[u,v]): zero(a-b)
relations = [x+y-u, lam*x-y-v,
             lam*(lam+1)*x*x+(lam+1)*y*y-lam*u*u-v*v]
for rel in relations: zero(rel.subs({u:G[0],v:G[1]}, simultaneous=True))
# Inverse output coordinates give the upper-bound differential equation.
C,D = s.Function('C'),s.Function('D')
mixed = s.diff(C(G[0])+D(G[1]), x, y)
expected = s.Subs(s.diff(C(t),t,2),t,G[0])-lam*s.Subs(s.diff(D(t),t,2),t,G[1])
zero(mixed-expected)
a,b,d = s.symbols('a b d', real=True)
Cg=lam*d*u*u/2+a*u
Dg=d*v*v/2+b*v
Ag=-(lam+1)*(lam*d*x*x/2+(a+lam*b)*x/(lam+1))
Bg=-(lam+1)*d*y*y/2-(a-b)*y
zero((Ag+Bg+Cg+Dg).subs({u:G[0],v:G[1]}, simultaneous=True))
# Three coefficients are independent modulo componentwise constants: derivatives
# of C determine (a,d), and then the derivative of D determines b.
coef=s.Matrix([[1,0,0],[0,0,lam],[0,1,0]])
assert s.factor(coef.det()) == -lam
print('FAMILY lambda>0: inverse, determinant -(lambda+1), three relations,')
print('  mixed equation C_second(u)-lambda*D_second(v)=0, and 3-parameter solution verified')
print('  independence coefficient determinant = -lambda (nonzero for lambda>0)')

# Chain-rule cancellation under exactly Remark 5's coordinate changes.
ux,uy,vx,vy,ap,bp,cp,dp = s.symbols('ux uy vx vy ap bp cp dp', nonzero=True)
zero(((cp*ux/ap)*(dp*vy/bp))/((cp*uy/bp)*(dp*vx/ap)) - ux*vy/(uy*vx))
zero(J[0,0]*J[1,1]/(J[0,1]*J[1,0])+1/lam)
print('INVARIANT: R_tilde(a(x),b(y))=R(x,y); Gaussian R=-1/lambda')

# Fresh computation for every member of the exact Koudou--Vallois list.
functions = [alpha/t,
             1/(alpha*(s.exp(beta*t)-1)),
             s.log((1+beta*t)/(beta*t))/alpha,
             s.log((s.exp(beta*t)+delta-1)/(s.exp(beta*t)-1))/alpha]
points = [(1,1,2),(s.log(2)/beta,s.log(2)/beta,s.log(4)/beta),
          (1/beta,1/beta,2/beta),(s.log(2)/beta,s.log(2)/beta,s.log(4)/beta)]
for i,(f,(xx,yy1,yy2)) in enumerate(zip(functions,points),1):
    fp=s.diff(f,t)
    R=fp.subs(t,x+y)/(fp.subs(t,x+y)-fp.subs(t,x))
    r1=s.factor(s.simplify(s.expand_log(R.subs({x:xx,y:yy1},simultaneous=True),force=True)))
    r2=s.factor(s.simplify(s.expand_log(R.subs({x:xx,y:yy2},simultaneous=True),force=True)))
    difference=s.factor(r1-r2)
    assert difference != 0
    numerator,denominator=s.fraction(difference)
    # All differences are strictly negative for positive parameters. Check signs
    # by exact polynomial coefficients, not sampling or floating point.
    assert all(c < 0 for c in s.Poly(numerator,alpha,beta,delta).coeffs())
    assert all(c > 0 for c in s.Poly(denominator,alpha,beta,delta).coeffs())
    print(f'KV family {i}: R at two positive-domain points = {r1}, {r2}')
    print(f'  difference = {difference} < 0 for alpha,beta,delta>0')

# General constant-invariant obstruction, cleared only where derivatives are nonzero.
r,p,q=s.symbols('r p q')
zero(s.expand(r*(q-p)-q)-(-r*p+(r-1)*q))
print('OBSTRUCTION: constant r gives r*p=(r-1)*q; r=1 forces p=0.')
print('  For r!=1, common s>x1,x2 forces f_prime(x1)=f_prime(x2).')
print('  Hence f is affine, incompatible with a decreasing bijection (0,infinity)->(0,infinity).')

# At lambda=1 the unnormalized standard Gaussian density transforms with |J|=2.
exponent=s.factor((x*x+y*y).subs(inv,simultaneous=True).subs(lam,1)/2)
zero(exponent-(u*u+v*v)/4)
print('GAUSSIAN lambda=1: exponent (u^2+v^2)/4 and absolute Jacobian 2')
assert not any(isinstance(a,s.Float) for obj in [J,exponent]+functions for a in obj.atoms())
print('ALL_CHECKS_OK')
