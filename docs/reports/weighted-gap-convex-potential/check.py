#!/usr/bin/env python3
"""Exact shooting certificates for Ahrami Open Problem 3.3.

Only Python integers, fractions and exact SymPy expressions are used.
All interval arithmetic rounds OUTWARD to a rational grid. See README.md
for the analytic remainder estimate and the eigenvalue-order bridge.
"""
from fractions import Fraction as Q
import sympy as sp

GRID=10**60

class I:
 def __init__(self,a,b=None): self.a=Q((Q(a)*GRID).__floor__(),GRID); self.b=Q((Q(a if b is None else b)*GRID).__ceil__(),GRID); assert self.a<=self.b
 def __add__(self,o):
  o=o if isinstance(o,I) else I(o);return I(self.a+o.a,self.b+o.b)
 __radd__=__add__
 def __neg__(self): return I(-self.b,-self.a)
 def __sub__(self,o): return self+-asI(o)
 def __rsub__(self,o): return asI(o)+-self
 def __mul__(self,o):
  o=asI(o);v=[self.a*o.a,self.a*o.b,self.b*o.a,self.b*o.b];return I(min(v),max(v))
 __rmul__=__mul__
 def __truediv__(self,q): assert Q(q)>0;return I(self.a/Q(q),self.b/Q(q))
 def absmax(self):return max(abs(self.a),abs(self.b))
 def __repr__(self):return f'[{self.a}, {self.b}]'
def asI(x):return x if isinstance(x,I) else I(x)

def atan(q,N):
 s=sum(((-1)**k*Q(1,(2*k+1)*q**(2*k+1)) for k in range(N)),Q(0))
 e=(-1)**N*Q(1,(2*N+1)*q**(2*N+1))
 return I(min(s,s+e),max(s,s+e))
PI=16*atan(5,90)-4*atan(239,30)
eps=Q(1,100);K=160
# For |t| <= 2: |A+B*t| <= 64. The sum of absolute Taylor
# coefficients weighted by 2^k is <= 2 exp(16) < 2*3^16.
TAIL=Q(2*3**16,2**(K+1))

def det(lam,potential):
 A=-(PI*PI)*lam*(1-eps*PI/2)
 B=eps*PI*PI*PI*((1 if potential else 0)-lam)
 assert A.absmax()+2*B.absmax()<64
 c0=I(0);c1=I(1);coeff=[c0,c1]
 for k in range(K-1):
  coeff.append((A*coeff[k]+(B*coeff[k-1] if k else 0))/((k+1)*(k+2)))
 s=sum(coeff,I(0));return s+I(-TAIL,TAIL)

def sign(lam,p):
 d=det(lam,p)
 assert d.a>0 or d.b<0,(lam,p,d)
 return 1 if d.a>0 else -1

def bracket(n,p):
 a=Q(9,10) if n==1 else Q(39,10)
 b=Q(11,10) if n==1 else Q(41,10)
 sa=sign(a,p);assert sign(b,p)==-sa
 for _ in range(36):
  m=(a+b)/2
  if sign(m,p)==sa:a=m
  else:b=m
 scale=10**12
 a=Q((a*scale).__floor__(),scale);b=Q((b*scale).__ceil__(),scale)
 assert sign(a,p)==sa and sign(b,p)==-sa
 return I(a,b)


def fidelity():
    x = sp.symbols('x', real=True)
    lam = sp.symbols('lambda', real=True)
    u = sp.Function('u')(x)
    G = sp.Function('G')(x)
    V = sp.Function('V')(x)
    w = sp.Function('w')(x)
    q = lam*w-V
    flux = G*(sp.diff(u,x)**2+q*u**2)-sp.diff(G,x)*u*sp.diff(u,x)+sp.diff(G,x,2)*u**2/2
    integrand = (2*sp.diff(G,x)*q+G*sp.diff(q,x)+sp.diff(G,x,3)/2)*u**2
    residual = sp.expand(sp.diff(flux,x)-integrand).subs(sp.diff(u,x,2),-q*u)
    assert sp.simplify(residual) == 0
    # Dirichlet endpoints make the boundary term G' u u' vanish.
    for n in [1,2]:
        un = sp.sqrt(2/sp.pi)*sp.sin(n*x)
        assert sp.simplify(-sp.diff(un,x,2)-n*n*un) == 0
        assert un.subs(x,0) == 0 and un.subs(x,sp.pi) == 0
        assert sp.integrate(un**2,(x,0,sp.pi)) == 1
        assert sp.simplify(sp.pi*sp.diff(un,x).subs(x,sp.pi)**2) == 2*n*n
    print('FIDELITY: Section 2 hypervirial identity and its normalized Dirichlet cases n=1,2: OK')
    S = {}
    for n in [1,2]:
        y = x-sp.pi/2
        sn = sp.sin(n*x)
        v = ((sp.pi/2)**2-y**2)*sp.cos(n*x)/(4*n)+y*sn/(4*n*n)
        assert sp.simplify(sp.diff(v,x,2)+n*n*v-y*sn) == 0
        assert sp.simplify(v.subs(x,0)) == 0 and sp.simplify(v.subs(x,sp.pi)) == 0
        assert sp.integrate(v*sn,(x,0,sp.pi)) == 0
        S[n] = sp.simplify(2/sp.pi*sp.integrate(y*sn*v,(x,0,sp.pi)))
        assert sp.simplify(S[n]-(n*n*sp.pi**2-15)/(48*n**4)) == 0
    assert sp.simplify(S[1]-7*S[2]+(45+4*sp.pi**2)/256) == 0
    print('PERTURBATION_COEFFICIENT: -(45+4*pi^2)/256: OK (not used for the fixed certificate)')

fidelity()
assert PI.b-PI.a < Q(1,10**58)
assert Q(3) < PI.a and PI.b < Q(22,7)
# Integer verification of the angle in Machin's identity.
z = sp.expand((5+sp.I)**4*(239-sp.I))
assert sp.re(z) == sp.im(z) > 0
assert Q(1,2) < Q(49,50) and Q(51,50) < Q(3,2)
sx = sp.Symbol('x')
assert sp.diff(sx/100,sx,2) == 0
assert sp.diff(1+(sx-sp.pi/2)/100,sx,2) == 0
assert Q(1)-eps*PI.b/2 > Q(49,50)
assert Q(1)+eps*PI.b/2 < Q(51,50)
assert eps*PI.b < Q(1,25) < 2
print('CLASS: M=2; Nlo=1/2; Nhi=3/2; V in [0,1/25]; w in (49/50,51/50); affine: OK')

print('PI_WIDTH_LE_10^-58',PI.b-PI.a<Q(1,10**58))
print('EPSILON',eps,'TAYLOR_DEGREE',K,'TAIL_LE_10^-39',TAIL<Q(1,10**39))
roots={}
for p in [False,True]:
 for n in [1,2]:
  roots[p,n]=bracket(n,p);print('ROOT', 'V=eps*x' if p else 'V=0',n,roots[p,n])
  for side,x in [('left',roots[p,n].a),('right',roots[p,n].b)]:
   d=det(x,p)
   sc=10**18
   dlo=Q((d.a*sc).__floor__(),sc);dhi=Q((d.b*sc).__ceil__(),sc)
   assert dlo>0 or dhi<0
   print('DETERMINANT',side,I(dlo,dhi))
# n^2/(51/50) <= lambda_n <= (n^2+1/25)/(49/50), using min-max.
assert Q(52,49)<Q(200,51)<Q(150,17)
for p in [False,True]:
 assert roots[p,1].a>0 and roots[p,1].b<Q(200,51)
 assert Q(52,49)<roots[p,2].a and roots[p,2].b<Q(150,17)
print('ORDER_LOCALIZATION: lambda_1 <= 52/49; lambda_2 >= 200/51; lambda_3 >= 150/17: OK')
diff=roots[True,2]-roots[True,1]-roots[False,2]+roots[False,1]
print('DIFF',diff)
assert diff.b<Q(-3,100000)
print('DIFFERENCE_UPPER_BOUND',Q(-3,100000))
print('ALL_CHECKS_OK')
