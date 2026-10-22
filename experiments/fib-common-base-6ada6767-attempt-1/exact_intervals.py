"""Outward fixed-grid rational intervals and forward derivatives.
Arithmetic follows unequal-diagonal M-pair Definition 4.1.
"""
from fractions import Fraction as Q
GRID=10**120
NVAR=1
def ceildiv(a,b):return -((-a)//b)
class I:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None):self.lo=lo;self.hi=lo if hi is None else hi;assert self.lo<=self.hi
    @classmethod
    def enc(cls,x):
        x=Q(x)*GRID;return cls(x.numerator//x.denominator,ceildiv(x.numerator,x.denominator))
    @classmethod
    def box(cls,c,rho):return cls(I.enc(c-rho).lo,I.enc(c+rho).hi)
    def __add__(self,y):
        y=asI(y);return I(self.lo+y.lo,self.hi+y.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,y):return self+-asI(y)
    def __rsub__(self,y):return asI(y)+-self
    def __mul__(self,y):
        y=asI(y);v=[a*b for a in [self.lo,self.hi] for b in [y.lo,y.hi]];return I(min(v)//GRID,ceildiv(max(v),GRID))
    __rmul__=__mul__
    def inv(self):
        assert self.lo*self.hi>0,(self.lo,self.hi)
        return I(GRID*GRID//self.hi,ceildiv(GRID*GRID,self.lo))
    def __truediv__(self,y):return self*asI(y).inv()
    def __rtruediv__(self,y):return asI(y)*self.inv()
    def mag(self):return Q(max(abs(self.lo),abs(self.hi)),GRID)
    def bounds(self):return [str(Q(self.lo,GRID)),str(Q(self.hi,GRID))]
    def rounded(self,n=12):
        g=10**n;return [str(Q((self.lo*g)//GRID,g)),str(Q(ceildiv(self.hi*g,GRID),g))]
    def nz(self):return self.lo*self.hi>0

def asI(y):return y if isinstance(y,I) else I.enc(y)
class AD:
    __slots__=('v','d')
    def __init__(self,v,d=None):self.v=asI(v);self.d=[I.enc(0) for k in range(NVAR)] if d is None else d
    def __add__(self,y):
        y=asAD(y);return AD(self.v+y.v,[a+b for a,b in zip(self.d,y.d)])
    __radd__=__add__
    def __neg__(self):return AD(-self.v,[-a for a in self.d])
    def __sub__(self,y):return self+-asAD(y)
    def __rsub__(self,y):return asAD(y)+-self
    def __mul__(self,y):
        y=asAD(y);return AD(self.v*y.v,[a*y.v+self.v*b for a,b in zip(self.d,y.d)])
    __rmul__=__mul__
    def inv(self):return AD(1/self.v,[-a/(self.v*self.v) for a in self.d])
    def __truediv__(self,y):return self*asAD(y).inv()
    def __rtruediv__(self,y):return asAD(y)*self.inv()
def asAD(y):return y if isinstance(y,AD) else AD(y)
def independent(box):return [AD(x,[I.enc(int(i==j)) for j in range(len(box))]) for i,x in enumerate(box)]
def contraction(fun,centres,C,rho):
    global NVAR
    NVAR=len(centres);box=[I.box(c,rho) for c in centres]
    f0=fun([I.enc(c) for c in centres]);df=fun(independent(box));n=NVAR
    cf=[sum((C[i][k]*f0[k] for k in range(n)),I.enc(0)) for i in range(n)]
    E=[[I.enc(int(i==j))-sum((C[i][k]*df[k].d[j] for k in range(n)),I.enc(0)) for j in range(n)] for i in range(n)]
    eta=max(a.mag() for a in cf);qs=[sum(e.mag() for e in row) for row in E];q=max(qs)
    assert q<1 and eta+q*rho<rho,(str(eta),str(q),str(rho))
    return box,{'eta':str(eta),'q_rows':list(map(str,qs)),'q':str(q),'rho':str(rho),'inclusion_margin':str(rho-eta-q*rho),'grid_denominator':str(GRID)}
