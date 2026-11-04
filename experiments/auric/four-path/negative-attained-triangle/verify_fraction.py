"""Standard-library-only exact rational verifier; no Arb import or floats.
Verifies the finite binary partition and its polynomial positivity tests.
"""
import sys,json,time,hashlib
from fractions import Fraction as F
from math import isqrt
assert sys.version_info >= (3,12)
class I:
    def __init__(self,l,h=None):
        if isinstance(l,I):self.l,self.h=l.l,l.h
        else:self.l,self.h=F(l),F(l if h is None else h)
        assert self.l<=self.h
    def __add__(self,o):
        o=ii(o);return I(self.l+o.l,self.h+o.h)
    __radd__=__add__
    def __neg__(self):return I(-self.h,-self.l)
    def __sub__(self,o):return self+-ii(o)
    def __rsub__(self,o):return ii(o)+-self
    def __mul__(self,o):
        o=ii(o);v=[self.l*o.l,self.l*o.h,self.h*o.l,self.h*o.h];return I(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=ii(o);assert o.l*o.h>0;return self*I(1/o.h,1/o.l)
    def __rtruediv__(self,o):return ii(o)/self
    def __pow__(self,n):
        assert isinstance(n,int) and n>=0
        if n==0:return I(1)
        if n%2:return I(self.l**n,self.h**n)
        return I(0 if self.l<=0<=self.h else min(self.l**n,self.h**n),max(self.l**n,self.h**n))
    def sqrt(self):
        assert self.l>=0;mesh=10**40
        def low(v):return F(isqrt(v.numerator*mesh*mesh//v.denominator),mesh)
        lo=low(self.l);hi=low(self.h)
        if hi*hi<self.h:hi+=F(1,mesh)
        assert lo*lo<=self.l and hi*hi>=self.h
        return I(lo,hi)
    def meet(self,o):return I(max(self.l,o.l),min(self.h,o.h))
    def positive(self):return self.l>0
    def serial(self):return [str(self.l),str(self.h)]
def ii(x):return x if isinstance(x,I) else I(x)
class D:
    def __init__(self,v,g=None):self.v=ii(v);self.g=g or (I(0),)*3
    def __add__(self,o):
        o=dd(o);return D(self.v+o.v,tuple(a+b for a,b in zip(self.g,o.g)))
    __radd__=__add__
    def __neg__(self):return D(-self.v,tuple(-a for a in self.g))
    def __sub__(self,o):return self+-dd(o)
    def __rsub__(self,o):return dd(o)+-self
    def __mul__(self,o):
        o=dd(o);return D(self.v*o.v,tuple(a*o.v+self.v*b for a,b in zip(self.g,o.g)))
    __rmul__=__mul__
    def __pow__(self,n):
        if n==0:return D(1)
        return D(self.v**n,tuple(n*self.v**(n-1)*a for a in self.g))
def dd(x):return x if isinstance(x,D) else D(x)
# I delegates binary operations involving D through NotImplemented.
for opname in ['__add__','__sub__','__mul__']:
    old=getattr(I,opname)
    def delegated(self,other,method=old):
        if isinstance(other,D):return NotImplemented
        return method(self,other)
    setattr(I,opname,delegated)
I.__radd__=I.__add__;I.__rmul__=I.__mul__
h=I(F(1,2)).sqrt();s=((1-h)/2).sqrt();c=((1+h)/2).sqrt()
e2=4*(h+s)**2;d2=4*(h+c)**2;k=h-s

def polynomials(t,a,b,mode):
    u=a;v=(1-a)*b
    r=1-t*u;z=1-t*v
    Hcoef=I(F(4,25))+I(F(2,3))*t
    HH=t*Hcoef
    delta=I(F(37,150))+I(F(2,3))*t if mode==0 else I(F(17,60))+I(F(1,3))*t
    A0=e2-HH;B0=d2-HH-delta
    rho=4*s*h-2*s*k*u+2*h*k*v-t*(s*h-s*k*u*u+h*k*v*v)
    W=4*rho-Hcoef
    Cbase=c*c*u*(2-t*u)+h*h*v*(2-t*v)
    P=W*B0-4*A0*Cbase
    N=(1-a)**2*a*b*(1-b)*(2-t*u)*(2-t*v)*(4-t*(1+u+v))*((2-t)**2-t*t*(u-v)**2)
    Q=r*r*z*z*P*P-16*A0*A0*c*c*h*h*t*N
    return W,B0,P,Q

def check(box,mode):
    values=[I(l,r) for l,r in box];raw=polynomials(*values,mode)
    if all(z.positive() for z in raw):return [z.l for z in raw]
    dual=[D(x,tuple(I(int(i==j)) for j in range(3))) for i,x in enumerate(values)]
    dv=polynomials(*dual,mode);cv=polynomials(*[I((l+r)/2) for l,r in box],mode)
    out=[]
    for direct,d,center in zip(raw,dv,cv):
        mean=center
        for derivative,(l,r) in zip(d.g,box):mean+=derivative*I(-(r-l)/2,(r-l)/2)
        bound=direct.meet(mean);assert bound.positive(),{'mode':mode,'box':box,'bound':bound.serial()};out.append(bound.l)
    return out

def auxiliary():
    Q=I(F(17,500));R=s*I(F(407,1000))+h
    C=((c*c+h*h)*(2-I(F(139,100))**2/2)).sqrt()
    AA=e2-4*R*R;BB=d2-4*C*C
    low=AA+(BB-AA)*Q-8*(Q*(1-Q)).sqrt()*R*C
    assert R.l>F(4,5) and C.l>1 and (BB-AA).h<4
    assert low.l>F(2,5)
    Hmax=I(F(61,25))/25+I(F(61,25))**2/24
    assert Hmax.h<F(2,5)
    # Sealed upper-comparator argument: exact quadratic endpoint tests.
    D0=4*(2*h*s+s*s-h*c)
    P0=I(F(8,25))-2*D0
    P1=I(F(18,25))+I(F(8,3))-3*D0-4*h*c
    P2=I(F(84,25))-6*h*c
    assert P2.h<0 and P0.l>0 and (P0+P1*I(F(2,3))+P2*I(F(4,9))).l>0
    assert (e2-I(F(4,25))*I(F(593,1000))-I(F(2,3))*I(F(593,1000))**2).l>4
    # Independently reproduce the already known J-only failure.
    xi=I(F(59,25));q=I(F(1,5));p=1-q;x=I(F(51,100));y=(1-x*x).sqrt();L=1+x;R0=s*x+h
    mu=(s*k*p.sqrt()/R0+c*q.sqrt()/y)/2
    nu=(k*p.sqrt()+(c*c*x-h*h)*q.sqrt()/(c*y))/2
    tc=h*q.sqrt()*L.sqrt()/x
    omega=nu-tc*tc/(4*mu);assert omega.h<0
    gaincoef=p.sqrt()*k*k/(8*R0)+2*h*q.sqrt()*L/(3*I(3).sqrt()*x)
    J=omega**2/(4*gaincoef);gamma=p.sqrt()*R0+q.sqrt()*c*y
    can=e2*p+d2*q-4*gamma**2;loss=8*gamma*J+4*J**2;glow=can-loss
    expected={'omega':(-354173,-354172),'J':(83979,83980),'Gamma0':(1162419,1162420),'Gcan':(523676,523677),'loss':(809168,809169),'Glow':(-285492,-285491)}
    values={'omega':omega,'J':J,'Gamma0':gamma,'Gcan':can,'loss':loss,'Glow':glow}
    for name,z in values.items():
        l,r=expected[name];assert F(l,10**6)<=z.l<=z.h<=F(r,10**6),name
    assert glow.h<F(-1,5)
    return {'small_slab_lower':low.serial(),'Hmax':Hmax.serial(),'comparator_P0':P0.serial(),'comparator_P_at_2over3':(P0+P1*I(F(2,3))+P2*I(F(4,9))).serial(),'known_J_only_failure':{n:z.serial() for n,z in values.items()}}

def main():
    start=time.monotonic();cert=json.load(open(sys.argv[1]));assert cert['status']=='complete';trees=cert['trees'];assert len(trees)==2
    roots=[[(F(0),F(11,100)),(F(0),F(1)),(F(0),F(1))],[(F(11,100),F(593,1000)),(F(0),F(1)),(F(0),F(1))]]
    leaves=0;minima=[None]*4
    for mode,tree in enumerate(trees):
        stack=[(roots[mode],tree)]
        while stack:
            box,node=stack.pop()
            if node==1:
                lower=check(box,mode);leaves+=1
                for j,v in enumerate(lower):minima[j]=v if minima[j] is None else min(minima[j],v)
            else:
                assert isinstance(node,list) and len(node)==3
                j,left,right=node;assert type(j) is int and 0<=j<3
                l,r=box[j];m=(l+r)/2;lb=list(box);rb=list(box);lb[j]=(l,m);rb[j]=(m,r)
                stack.extend([(lb,left),(rb,right)])
    assert leaves==cert['leaves']
    clean_bounds=[F(12,5),F(49,5),F(1,10),F(3,500)]
    assert all(v>b for v,b in zip(minima,clean_bounds))
    out={'status':'verified','arithmetic':'exact fractions; square roots bracketed on rational mesh 1e-40 by integer isqrt','leaves':leaves,'polynomial_lower_bounds':[str(x) for x in minima],'simple_strict_bounds':dict(zip(['W','B0','P','Q'],map(str,clean_bounds))),'auxiliary':auxiliary()}
    with open(sys.argv[2],'w') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({'status':out['status'],'leaves':leaves,'seconds':time.monotonic()-start,'minimum_positive':all(x>0 for x in minima)}))
if __name__=='__main__':main()
