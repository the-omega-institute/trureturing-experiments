#!/usr/bin/env python3
"""Exact, outward-rounded rational checks. No floating point is used.

Analytic inputs: alternating-series bounds for atan/cos, the explicit Robin
quantization equation, separation of variables, and the lower rectangle rule
for a decreasing nonnegative function. These finite checks do not prove the
conjecture or a uniform asymptotic estimate.
"""
from fractions import Fraction as Q
from math import isqrt, factorial
from functools import lru_cache
import sympy as sp

S = 10**28
class I:
    __slots__ = ('a', 'b')
    def __init__(self, x=0, y=None):
        if isinstance(x,I): self.a,self.b=x.a,x.b; return
        x,y=Q(x),Q(x if y is None else y)
        self.a=x.numerator*S//x.denominator
        self.b=-((-y.numerator*S)//y.denominator)
        assert self.a<=self.b
    @classmethod
    def raw(cls,a,b):
        z=object.__new__(cls);z.a,z.b=a,b;assert a<=b;return z
    def __add__(self,x):
        x=I(x);return I.raw(self.a+x.a,self.b+x.b)
    __radd__=__add__
    def __neg__(self):return I.raw(-self.b,-self.a)
    def __sub__(self,x):return self+-I(x)
    def __rsub__(self,x):return I(x)+-self
    def __mul__(self,x):
        x=I(x);v=[self.a*x.a,self.a*x.b,self.b*x.a,self.b*x.b]
        return I.raw(min(v)//S,-((-max(v))//S))
    __rmul__=__mul__
    def __truediv__(self,x):
        x=I(x);assert not x.a<=0<=x.b
        v=[Q(u*S,w) for u in (self.a,self.b) for w in (x.a,x.b)]
        a,b=min(v),max(v)
        return I.raw(a.numerator//a.denominator,-((-b.numerator)//b.denominator))
    def __rtruediv__(self,x):return I(x)/self
    def __pow__(self,n):
        assert isinstance(n,int) and n>=0
        z=I(1)
        for _ in range(n):z=z*self
        return z
    def sqrt(self):
        assert self.a>=0
        a,b=isqrt(self.a*S),isqrt(self.b*S)
        return I.raw(a,b+(b*b<self.b*S))
    def positive_part(self):return I.raw(max(self.a,0),max(self.b,0))
    def mid(self):return Q(self.a+self.b,2*S)
    def absmax(self):return Q(max(abs(self.a),abs(self.b)),S)
    def show(self,digits=10):
        # Decimal endpoints rounded OUTWARDS by integer arithmetic.
        p=10**digits; a=self.a*p//S;b=-((-self.b*p)//S)
        def f(n):return ('-' if n<0 else '')+str(abs(n)//p)+'.'+str(abs(n)%p).zfill(digits)
        return '['+f(a)+', '+f(b)+']'

@lru_cache(None)
def atan_small(q):
    assert 0<=q<=Q(1,2)
    term=q;z=Q(0)
    for k in range(56):
        z+=(-1)**k*term/(2*k+1);term*=q*q
    r=term/113
    return I(z-r,z+r)
PI = 16*atan_small(Q(1,5))-4*atan_small(Q(1,239))

def atan_point(q):
    q=Q(q)
    if q<0:return -atan_point(-q)
    if q>1:return PI/2-atan_point(1/q)
    if q>Q(1,2):return atan_small(Q(1,2))+atan_small((q-Q(1,2))/(1+q/2))
    return atan_small(q)
def atan(x):
    x=I(x);a=atan_point(Q(x.a,S));b=atan_point(Q(x.b,S))
    return I.raw(a.a,b.b)

def cos(x):
    x=I(x)
    n=(x.mid()/PI.mid()/2+Q(1,2)).numerator//(x.mid()/PI.mid()/2+Q(1,2)).denominator
    x=x-2*n*PI;assert x.absmax()<4
    term=I(1);z=I(1)
    for k in range(1,70):term=term*(-x*x)/((2*k-1)*(2*k));z=z+term
    rem=Q(4**140,factorial(140))
    return z+I(-rem,rem)

def powhalf(x,a):
    a=Q(a);assert a>=0 and a.denominator in (1,2)
    x=I(x);return (x**(a.numerator//a.denominator))*(x.sqrt() if a.denominator==2 else 1)

def gammahalf(x):
    x=Q(x);assert x>0 and x.denominator in (1,2)
    if x.denominator==1:return I(factorial(x.numerator-1))
    k=(x.numerator-1)//2
    return I(Q(factorial(2*k),4**k*factorial(k)))*PI.sqrt()

def L(a,d):return gammahalf(a+1)/(powhalf(4*PI,Q(d,2))*gammahalf(a+1+Q(d,2)))

def boundary(a,b):
    a=Q(a);b=I(b);p=1+b*b
    if a.denominator==1:
        integral=atan(1/b);moment=I(1);start=Q(0)
    else:
        integral=PI/2*(p.sqrt()-b);moment=PI/4;start=Q(1,2)
    while start<a:
        integral=p*integral-b*moment
        start+=1;moment=moment*(2*start)/(2*start+1)
    return 4*integral/PI-1

def beta_bracket(a):
    lo,hi=Q(1,100),Q(2)
    assert boundary(a,lo).a>0 and boundary(a,hi).b<0
    for _ in range(42):
        m=(lo+hi)/2;r=boundary(a,m)
        if r.a>0:lo=m
        elif r.b<0:hi=m
        else:raise AssertionError('insufficient rational precision')
    return I(lo,hi)

def frequencies(t,b,l=1):
    t,b,l=I(t),I(b),Q(l);eta=b*t
    ans=[];j=1
    while (j-1)*Q(PI.a,S)<l*Q(t.b,S):
        lo=Q(0) if j==1 else (j-1)*Q(PI.a,S)/l
        hi=j*Q(PI.b,S)/l
        for _ in range(46):
            m=(lo+hi)/2
            f=l*I(m)+2*atan(I(m)/eta)-j*PI
            if f.b<0:lo=m
            elif f.a>0:hi=m
            else:break
        assert hi-lo<Q(1,10**11)
        ans.append(I(lo,hi));j+=1
    return ans

def riesz(a,t,b):
    t=I(t)
    return sum((powhalf((t*t-k*k).positive_part(),a) for k in frequencies(t,b)),I(0))

def periodic(a,theta,M=768):
    # Closed Bernoulli formula for integer exponents 1 and 2.
    u=theta/PI
    if a==1:return -PI**2*(u*u-u+Q(1,6))
    if a==2:return -2*PI**3/3*(u**3-Q(3,2)*u*u+u/2)
    ans=I(0)
    for m in range(1,M+1):
        ans=ans+cos(2*m*theta-PI*(a+1)/2)/powhalf(m,a+1)
    tail=1/(I(a)*powhalf(M,a))
    return ans+I(-Q(tail.b,S),Q(tail.b,S))

def fidelity():
    # Lemma 7.2(ii): beta_W(1/2,0)=3/4.
    # Rational half-angle substitution s=2u/(1+u^2) in its integral.
    u=sp.symbols('u',real=True);b=sp.Rational(3,4)
    integrand=2*b*(1-u*u)**2/((1+u*u)*(b*b*(1+u*u)**2+4*u*u))
    anti=-2*b*sp.atan(u)+sp.Rational(5,4)*(sp.atan(3*u)+sp.atan(u/3))
    assert sp.cancel(sp.diff(anti,u)-integrand)==0
    # atan(3)+atan(1/3)=pi/2, atan(1)=pi/4: integral=pi/4.
    assert -2*b*sp.Rational(1,4)+sp.Rational(5,4)*sp.Rational(1,2)==sp.Rational(1,4)
    print('FIDELITY beta_W(1/2,0)=3/4: exact rational antiderivative checked')
    assert boundary(Q(1,2),Q(3,4)).a<=0<=boundary(Q(1,2),Q(3,4)).b
    count=0
    for eta in (1,3,10):
        # Here bt=eta, t=30; intervals contain all first eight roots.
        for j,k in enumerate(frequencies(30,Q(eta,30))[:8],1):
            lower=j*PI-2*atan(j*PI/eta)
            upper=j*PI-2*atan((j-1)*PI/eta)
            assert lower.b<k.a and k.b<upper.a
            count+=1
    print('FIDELITY Lemma 2.3: strict two-sided frequency bounds,',count,'cases')


def asymptotics():
    print('ASYMPTOTIC CHECKS b=3/4, theta=3*pi/4; intervals are rigorous finite readings')
    print('d=2 uses alpha=gamma+1/2 in the partial semiclassical limit, not a square-domain expansion')
    theta=3*PI/4;b=Q(3,4)
    for d in (1,2):
        for g in (Q(1),Q(2)):
            a=g+Q(d-1,2);G=gammahalf(a+1)/PI*periodic(a,theta)
            assert G.a>0
            print('PROFILE d=%d gamma=%s alpha=%s G=%s'%(d,g,a,G.show()))
            for n in (4,8,16):
                t=n*PI+theta-2*atan(1/I(b))
                R=riesz(a,t,b)
                norm=(R-L(a,1)*powhalf(t,2*a+1)-boundary(a,b)/2*powhalf(t,2*a))/powhalf(t,a)
                error=t*(norm-G)
                assert error.absmax()<10
                print('  n=%d t=%s normalized=%s t*(normalized-G)=%s'%(n,t.show(6),norm.show(),error.show(6)))
                if d==2:
                    ell=10**8
                    transverse_error=riesz(g,t,b)/(ell*L(g,1)*powhalf(t,a))
                    width=Q(transverse_error.b,S)
                    slab_norm=norm+I(-width,width)
                    print('    actual finite rectangle ell=%d normalized enclosure=%s'%(ell,slab_norm.show()))


def settlement_samples():
    print('STRICT THRESHOLD SAMPLES (finite checks; universal conjecture remains a research target)')
    for d in (1,2):
        for g in (Q(1),Q(2)):
            a=g+Q(d-1,2);bw=beta_bracket(a)
            # A rational parameter strictly ABOVE the certified Weyl-zero bracket.
            b=Q(bw.b,S)+Q(1,10**6)
            assert boundary(a,b).b<0
            t=4*PI+3*PI/4-2*atan(1/I(b));Ra=riesz(a,t,b)
            gap=Ra-L(a,1)*powhalf(t,2*a+1)
            assert gap.a>0
            if d==1:
                print('  d=1 gamma=%s betaW=%s b=%s t=%s R-W=%s'%(g,bw.show(),I(b).show(),t.show(6),gap.show()))
            else:
                # For each longitudinal mode f(p)=(t^2-k^2-p^2)_+^gamma,
                # sum_{j>=1} f(pi*j/ell) >= ell/pi*integral_0^inf f - f(0).
                # Transverse Robin eigenvalues are below (pi*j/ell)^2.
                # Thus R_rectangle/ell >= L(g,1)*R_(g+1/2) - R_g/ell.
                Rg=riesz(g,t,b);leading=L(g,2)*powhalf(t,2*g+2)
                ell=10**8
                lower=L(g,1)*Ra-Rg/ell
                certified_gap=lower-leading
                assert certified_gap.a>0
                const=L(g,1)*L(a,1)-L(g,2)
                assert const.a<=0<=const.b
                print('  d=2 gamma=%s betaW=%s b=%s t=%s ell=%d (R_rectangle/ell-W)_lower=%s'%(g,bw.show(),I(b).show(),t.show(6),ell,certified_gap.show()))
                print('    partial-limit error upper R_gamma/ell=%s'%(Rg/ell).show())

if __name__=='__main__':
    fidelity();asymptotics();settlement_samples();print('ALL_CHECKS_OK')
