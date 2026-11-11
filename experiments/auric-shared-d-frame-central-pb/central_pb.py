"""Exact full-M4 central-residual PB witnesses for degenerate elementary legs."""
import json
import sympy as s
r=s.symbols('r',real=True,nonzero=True)
J=lambda a,b:(a*b+b*a)/2
results=[]
for a,j in [(1,0),(2,1)]:
 for sigma in [1,-1]:
  ii=sigma*s.I
  K=lambda a,b:ii*(a*b-b*a)/2
  c=s.diag(1,2,3,0);t=s.zeros(4);t[j,3]=t[3,j]=1
  lam=a*(a*a-4*r*r)/(4*r*(a*a+2));x=s.zeros(4);x[j,j]=lam
  g=K(c,t)+r*s.eye(4)
  xi=s.zeros(4,1);xi[j]=r*(1-2*ii/a);xi[3]=1
  h=K(c,x);rr=J(c,x);kt=J(c,t)
  L=-2*J(x,h-kt)-2*J(t,rr+g)
  ell=s.factor((xi.conjugate().T*L*xi)[0]);Bs=[]
  for de in [1,-1]:
   R=x+de*(h-kt);S=t+de*(rr+g)
   rv=R*xi;sv=S*xi
   A=s.expand((rv.conjugate().T*rv)[0]);b=s.re(s.expand((rv.conjugate().T*sv)[0]));C=s.expand((sv.conjugate().T*sv)[0])
   Bs.append(s.factor(A*C-b*b))
  if a==1:
   wantB=(80*r**6+12*r**4+6*r**2+1)/576
   wantell=-(44*r**2+1)/12
   wantmargin=r**2*(41+962*r**2-40*r**4)/72
  else:
   wantB=(2*r**6+6*r**4+9*r**2+1)/9
   wantell=-2*(5*r**2+1)/3
   wantmargin=4*r**2*(1+19*r**2-2*r**4)/9
  assert all(s.simplify(B-wantB)==0 for B in Bs)
  assert s.simplify(ell-wantell)==0
  assert s.simplify(ell*ell-4*Bs[0]-wantmargin)==0
  # Large residual witness: x=0, normalized +1 eigenvector of t.
  xx=s.zeros(4);zz=s.zeros(4,1);zz[j]=zz[3]=1/s.sqrt(2)
  LL=-2*J(t,g);largeBs=[]
  for de in [1,-1]:
   R=-de*kt;S=t+de*g;rv=R*zz;sv=S*zz
   A=s.expand((rv.conjugate().T*rv)[0]);b=s.re(s.expand((rv.conjugate().T*sv)[0]));C=s.expand((sv.conjugate().T*sv)[0])
   largeBs.append(s.factor(A*C-b*b))
  assert all(s.simplify(B-s.Rational(a**4,16))==0 for B in largeBs)
  assert s.simplify((zz.conjugate().T*LL*zz)[0]+2*r)==0
  results.append({'a':a,'sigma':sigma,'small_B':str(wantB),'small_L':str(wantell),'small_squared_margin':str(wantmargin),'small_range':'0 < abs(r) <= a^2/4','large_B':str(s.Rational(a**4,16)),'large_L':'-2*r','large_range':'abs(r) > a^2/4'})
print(json.dumps({'full_M4_symbolic_checks':'passed','results':results},sort_keys=True))
