#!/usr/bin/env python3
"""Exact checks for the Duarte--da Mota--Skea rational-gradient candidate.

No bounded ansatz is used to assert universal nonexistence. The exact identities
below check the ingredients of the general paper argument stated in README.md.
Requires Python 3 and SymPy. All asserted arithmetic is symbolic or integral.
"""
import sympy as S
from itertools import product


def zero(label, expression):
    n = S.fraction(S.cancel(S.together(S.expand(expression))))[0]
    assert S.expand(n) == 0, (label, S.factor(n))
    print(label + ": OK")


x, y, z = S.symbols("x y z")
# Fidelity: source Section 3.3.1, harmonic oscillator y'' = -y.
# The source explicitly gives the conserved energy C1 = y**2 + z**2.
# The additional phase invariant is only a cross-check, not attributed to the source.
I1 = x - S.atan(y/z)
I2 = y*y + z*z
Dho = lambda f: S.diff(f,x) + z*S.diff(f,y) - y*S.diff(f,z)
zero("harmonic oscillator phase cross-check", Dho(I1))
zero("source harmonic oscillator energy C1", Dho(I2))
Snull=y/z
R=z
zero("source harmonic oscillator null form",Dho(Snull)-(1+Snull**2))
zero("source harmonic oscillator integrating factor",Dho(R)+R*Snull)
zero("source harmonic oscillator compatibility",S.diff(R,y)-S.diff(R,z)*Snull-S.diff(Snull,z)*R)

u,v,w,alpha,beta = S.symbols("u v w alpha beta", nonzero=True)
E = lambda f: u*S.diff(f,u) + alpha*v*S.diff(f,v) + beta*w*S.diff(f,w)
s=(u-1)/(u+1)
xx=s*s
AA=alpha*(v-1)/(v+1)/s
BB=beta*(w-1)/(w+1)/s
scale=s*(1-xx)
zero("normalization Dx = 1", E(xx)-scale)
a=lambda X: -1/(2*(1-X))
b=lambda X: -1/(2*X)
c=lambda h,X: h/(2*X*(1-X))
zero("Riccati A", E(AA)-scale*(a(xx)*AA**2+b(xx)*AA+c(alpha**2,xx)))
zero("Riccati B", E(BB)-scale*(a(xx)*BB**2+b(xx)*BB+c(beta**2,xx)))
for name, f in [("x",xx),("A",AA),("B",BB)]:
    zero("inversion invariance "+name, f.xreplace({u:1/u,v:1/v,w:1/w})-f)
zero("inversion reverses normalization", scale.xreplace({u:1/u})+scale)

J1=S.log(v)-alpha*S.log(u)
J2=S.log(w)-beta*S.log(u)
zero("first integral 1", E(J1))
zero("first integral 2", E(J2))
wedge=S.diff(J1,v)*S.diff(J2,w)-S.diff(J1,w)*S.diff(J2,v)
zero("independence dv wedge dw coefficient = 1/(vw)",wedge-1/(v*w))
assert wedge != 0

A,B,p,q=S.symbols("A B p q")
fA=a(x)*A*A+b(x)*A+c(p,x)
fB=a(x)*B*B+b(x)*B+c(q,x)
DD=lambda f:S.diff(f,x)+fA*S.diff(f,A)+fB*S.diff(f,B)
yy=A-B
zz=DD(yy)
zero("jet Jacobian determinant = 2ay",S.det(S.Matrix([yy,zz]).jacobian([A,B]))-2*a(x)*yy)
H=(z-b(x)*y-c(p,x)+c(q,x))/(a(x)*y)
As=(H+y)/2
Bs=(H-y)/2
subs={A:As,B:Bs}
zero("jet inverse y",yy.subs(subs)-y)
zero("jet inverse z",zz.subs(subs)-z)
Phi=3*z**2/(2*y)+((p-q)/(x*y)-1)*z/(x-1)+y**3/(8*(x-1)**2)+(3-2*(p+q+1)*x-x**2)*y/(8*x**2*(x-1)**2)+(p-q)/(2*x**2*(x-1))+(p-q)**2/(8*x**2*y*(x-1)**2)
zero("scalar jet Phi for symbolic p,q",DD(zz).subs(subs)-Phi)

# General elementary solution, checked through its two independent Riccati coordinates.
ss,U,C1,C2,gamma=S.symbols("s U C1 C2 gamma", nonzero=True)
V=S.symbols("V", nonzero=True)
ric=gamma/ss*(V-1)/(V+1)
ds=1/(2*ss)
dV=gamma*V/(ss*(1-ss**2))
zero("elementary solution Riccati",ds*S.diff(ric,ss)+dV*S.diff(ric,V)-(a(ss**2)*ric**2+b(ss**2)*ric+c(gamma**2,ss**2)))
zero("elementary U logarithmic derivative",S.diff((1+ss)/(1-ss),ss)*ds/((1+ss)/(1-ss))-1/(ss*(1-ss**2)))
assert S.diff(ric,V) != 0
print("two independent solution parameters: each Riccati coordinate depends nontrivially on its own parameter: OK")

# No-rational-form mechanism: for a closed horizontal eta=sum h_i dU_i,
# d(i_E eta)=0 gives E(h_j)+lambda_j h_j=0, hence E(U_j h_j)=0.
# This checks the exact differential identity, with closure explicitly imposed.
coords=(u,v,w)
weights=(S.Integer(1),alpha,beta)
h=[S.Function("h"+str(i))(*coords) for i in range(3)]
horizontal=sum(weights[i]*coords[i]*h[i] for i in range(3))
for j in range(3):
    residual=S.diff(horizontal,coords[j])-(E(h[j])+weights[j]*h[j])
    expected=sum(weights[i]*coords[i]*(S.diff(h[i],coords[j])-S.diff(h[j],coords[i])) for i in range(3))
    zero("Cartan closed-horizontal identity "+str(j),residual-expected)
    zero("constant coordinate coefficient identity "+str(j),E(coords[j]*h[j])-coords[j]*(E(h[j])+weights[j]*h[j]))
    inv=1/coords[j]
    zero("inversion negates dlog "+str(j),S.diff(inv,coords[j])/inv+1/coords[j])

# Exact Q-linear independence of 1,sqrt(2),m sqrt(3) for all m>=1:
# 1,sqrt(2),sqrt(3),sqrt(6) is a degree-four number-field basis;
# scaling one basis vector by a nonzero integer preserves independence.
t=S.Symbol("t")
minpoly=S.minimal_polynomial(S.sqrt(2)+S.sqrt(3),t)
assert S.Poly(minpoly,t).degree()==4
assert S.Poly(minpoly,t).is_irreducible
assert S.Poly(minpoly,t)==S.Poly(t**4-10*t**2+1,t)
print("nonresonant number-field basis: irreducible t^4-10t^2+1, degree 4: OK")
# Diagnostic finite check only: distinct monomial weights for m=1,...,5, degrees <=4.
exponents=list(product(range(5),repeat=3))
for m in range(1,6):
    weights_seen=set()
    for i,j,k in exponents:
        key=(i,j,m*k)
        assert key not in weights_seen
        weights_seen.add(key)
    zero("Phi member m="+str(m), (DD(zz).subs(subs)-Phi).subs({p:2,q:3*m*m}))
print("diagnostic weights: 125 distinct monomials per m=1,...,5: OK")

# A concrete generic point witnesses a nonempty denominator-safe chart at m=1.
s0=S.Rational(1,2);A0=S.Integer(2);B0=S.Integer(1)
x0=s0*s0;y0=A0-B0
z0=zz.subs({x:x0,A:A0,B:B0,p:2,q:3})
value=S.cancel(Phi.subs({x:x0,y:y0,z:z0,p:2,q:3}))
assert value.is_Rational and x0 not in (0,1) and y0 != 0
assert 2-x0*A0*A0 != 0 and 3-x0*B0*B0 != 0
print("denominator-safe m=1 point (x,y,z,Phi) =",(x0,y0,z0,value))
print("universal exclusion uses the coprime-eigenpolynomial and inversion argument in README; no finite search is its proof")
print("ALL_CHECKS_OK")
