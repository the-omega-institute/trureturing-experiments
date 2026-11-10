"""Exact algebra and rational interval checks for arXiv:2609.01373v1.

Uses only integers, fractions and SymPy exact expressions. The interval layer
rounds outwards to a rational grid. This checks formulas, not the conjecture.
"""
from fractions import Fraction as Q
from math import factorial, isqrt
import sympy as S

l, t, d = S.symbols('lambda tau D', real=True)
c = (1+l*l)/(1-l*l)
s = 2*l/(1-l*l)
a, b = 1-t+t*c, t*s
v = t+(1-t)*c
delta = S.factor(a*a-b*b)

def zero(x):
    assert S.cancel(x) == 0, S.factor(x)

def modulo_delta(x):
    num, _ = S.fraction(S.cancel(x))
    r = S.Poly(num, d).rem(S.Poly(d*d-delta, d)).as_expr()
    zero(r)

zero(delta - (1-l*l*(2*t-1)**2)/(1-l*l))
zero(a*c-b*s-v)
zero(v*v-((1-t)*s)**2-delta)
gamma = (1-d)/(1+d)
source_f = (1-gamma**2)*(1-l*l)*d/(1-(l*t-gamma*l*(1-t))**2)
modulo_delta(source_f - 2*d/(v+d))
kappa = b/(a+d)
modulo_delta((1+kappa*kappa)/(1-kappa*kappa)-a/d)
modulo_delta(2*kappa/(1-kappa*kappa)-b/d)
modulo_delta((1-l*l)*(1-kappa*kappa)/(1-l*kappa)**2 - 2*d/(v+d))
print('SYMBOLIC_OK: conditional determinant, Schmidt parameter, source squeezing fidelity, overlap identity')

I = S.eye(2)
X = S.Matrix([[0,1],[1,0]])
case_count = 0
for lam in [S.Rational(1,10), S.Rational(1,5), S.Rational(1,2), S.Rational(9,10)]:
    cc, ss = c.subs(l,lam), s.subs(l,lam)
    for tau in [S.Rational(0),S.Rational(1,4),S.Rational(1,2),S.Rational(3,4),S.Rational(1)]:
        w = S.sqrt(tau*(1-tau))
        for sign in [-1,1]:
            G = cc*I-sign*ss*X
            A = tau*G+(1-tau)*I
            E = w*(I-G)
            H = (1-tau)*G+tau*I
            V = ((1-tau)*G.inv()+tau*I)/2
            for x in H-E.T*A.inv()*E - V.inv()/2:
                zero(S.simplify(x))
            D = S.sqrt(A.det())
            assert D.is_positive
            for z in [S.Matrix([0,0]),S.Matrix([1,-2]),S.Matrix([2,3])]:
                m = -A.inv()*E*z
                cov_mean = (w*(I-G.inv())/2)*V.inv()*z
                for x in m-cov_mean: zero(S.simplify(x))
                # Hessian and completed center determine the entire conditional Gaussian.
                y0,y1 = S.symbols('y0 y1',real=True)
                y = S.Matrix([y0,y1])
                original = (y.T*A*y)[0]+2*(y.T*E*z)[0]+(z.T*H*z)[0]
                completed = ((y-m).T*A*(y-m))[0]+(z.T*(H-E.T*A.inv()*E)*z)[0]
                zero(S.expand(original-completed))
                centered_precision = A/D
                kk = S.simplify(sign*A[0,1]*(-1)/(A[0,0]+D))
                for x in centered_precision-((1+kk*kk)*I-sign*2*kk*X)/(1-kk*kk):
                    zero(S.simplify(x))
                # Actual pure Gaussian covariance fidelity, with canceled mean.
                target = S.diag(1,1,1,1)
                # Reordered coordinates (x1,x2,p1,p2), hence block diagonal.
                target = S.diag(G.inv()/2,G/2)
                output = S.diag(D*A.inv()/2,A/(2*D))
                gaussian_f = 1/S.sqrt((target+output).det())
                f_schmidt = (1-lam*lam)*(1-kk*kk)/(1-lam*kk)**2
                zero(S.simplify(gaussian_f-f_schmidt))
                gg = S.simplify((1-D)/(1+D))
                f_source = (1-gg*gg)*S.sqrt((1-lam*lam)*(1-lam*lam*(2*tau-1)**2))/(1-(lam*tau-gg*lam*(1-tau))**2)
                zero(S.simplify(gaussian_f-f_source))
                case_count += 1
print(f'EXACT_GRID_OK: {case_count} conditional states; lambda=1/10,1/5,1/2,9/10; five transmissions; both signs; three outcomes')

SCALE = 10**100

def floorq(q): return Q(q.numerator*SCALE//q.denominator,SCALE)
def ceilq(q): return -floorq(-q)
class B:
    def __init__(self,lo,hi=None):
        self.lo,self.hi=Q(lo),Q(lo if hi is None else hi)
        assert self.lo<=self.hi
    @staticmethod
    def of(x): return x if isinstance(x,B) else B(x)
    @staticmethod
    def out(lo,hi): return B(floorq(lo),ceilq(hi))
    def __add__(self,x):
        x=B.of(x);return B.out(self.lo+x.lo,self.hi+x.hi)
    __radd__=__add__
    def __neg__(self): return B(-self.hi,-self.lo)
    def __sub__(self,x): return self+-B.of(x)
    def __rsub__(self,x): return B.of(x)+-self
    def __mul__(self,x):
        x=B.of(x);q=[self.lo*x.lo,self.lo*x.hi,self.hi*x.lo,self.hi*x.hi];return B.out(min(q),max(q))
    __rmul__=__mul__
    def __truediv__(self,x):
        x=B.of(x);assert not x.lo<=0<=x.hi;return self*B.out(1/x.hi,1/x.lo)
    def __rtruediv__(self,x):return B.of(x)/self
    def __pow__(self,n):
        assert n>=0;r=B(1)
        for _ in range(n):r=r*self
        return r
    def contains(self,x):return self.lo<=x<=self.hi
    def shown(self,places=12):
        scale=10**places
        lo=self.lo.numerator*scale//self.lo.denominator
        hi=-((-self.hi.numerator*scale)//self.hi.denominator)
        def fmt(n):return ('-' if n<0 else '')+str(abs(n)//scale)+'.'+str(abs(n)%scale).zfill(places)
        return '['+fmt(lo)+', '+fmt(hi)+']'

def atan_small(x,n=160):
    assert 0<=x<=Q(1,2)
    terms=[(-1)**k*x**(2*k+1)/Q(2*k+1) for k in range(n)]
    a0=sum(terms,Q(0));a1=a0+(-1)**n*x**(2*n+1)/Q(2*n+1)
    return B.out(min(a0,a1),max(a0,a1))
PI=16*atan_small(Q(1,5))-4*atan_small(Q(1,239))
assert PI.lo>Q(3141592653589793238462643383279502884,10**36)
assert PI.hi<Q(3141592653589793238462643383279502885,10**36)

def trig(x,cos=False):
    assert max(abs(x.lo),abs(x.hi))<4
    n=48
    term=B(1) if cos else x
    total=term
    for k in range(1,n):
        denom=(2*k-1)*(2*k) if cos else (2*k)*(2*k+1)
        term=-term*x*x/denom;total=total+term
    power=2*n if cos else 2*n+1
    err=Q(4**power,factorial(power))
    return total+B(-err,err)

def qfun(p):
    x=PI*B.of(p);sn=trig(x);cs=trig(x,True)
    return sn*(1+sn)+x*cs

def hfun(p):
    sn=trig(PI*B.of(p));return 2*B.of(p)*sn/(1+sn)

lo,hi=Q(7,10),Q(3,4)
assert qfun(lo).lo>0 and qfun(hi).hi<0
for _ in range(64):
    mid=(lo+hi)/2;val=qfun(mid)
    if val.lo>0:lo=mid
    else:
        assert val.hi<0
        hi=mid
pstar=B(lo,hi)
score=hfun(pstar)
assert Q('0.7138035')<pstar.lo<pstar.hi<Q('0.7138045')
assert Q('0.6268285')<score.lo<score.hi<Q('0.6268295')
print('CERTIFIED_P_STAR',pstar.shown(18))
print('SOURCE_SCORE_0.626829',score.shown(18))
# Weak-squeezing exact source formula at the fully reflecting boundary.
for lam in [Q(1,10),Q(1,5)]:
    p=B(Q(1,2))+2*atan_small(lam)/PI
    val=(1-lam*lam)*p
    assert p.hi<pstar.lo
    assert (hfun(p)-val).contains(0)
    assert (hfun(p)-val).hi-(hfun(p)-val).lo<Q(1,10**55)
    print('SOURCE_WEAK_LAMBDA',str(lam),'P',p.shown(18),'OPTIMAL_SCORE',val.shown(18))
# Three strongly squeezed examples lie past the same certified stationary point.
# atan(lambda)>atan(1/2) for each lambda>1/2, and for 1/2 use its exact series.
p_half=B(Q(1,2))+2*atan_small(Q(1,2))/PI
assert p_half.lo>pstar.hi
for lam in [Q(1,2),Q(3,4),Q(9,10)]:
    print('SOURCE_PLATEAU_LAMBDA',str(lam),'OPTIMAL_SCORE',score.shown(18))
# q'(p)=2*pi*cos(pi*p)*(1+sin(pi*p))-pi^2*p*sin(pi*p),
# so it is strictly negative for 1/2<p<1; this is an analytic sign argument.
print('ALL_CHECKS_OK')
