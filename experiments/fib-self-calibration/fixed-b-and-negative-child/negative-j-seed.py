from fractions import Fraction as Fr
H=10**60

def ceildiv(a,b): return -((-a)//b)
class I:
    def __init__(self,lo,hi=None):
        self.lo=lo;self.hi=lo if hi is None else hi
        assert self.lo<=self.hi
    @staticmethod
    def exact(x):
        x=Fr(x)*H
        return I(x.numerator//x.denominator,ceildiv(x.numerator,x.denominator))
    def __add__(self,y):
        y=asI(y);return I(self.lo+y.lo,self.hi+y.hi)
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self,y): return self+-asI(y)
    def __rsub__(self,y): return asI(y)+-self
    def __mul__(self,y):
        y=asI(y);a=[x*z for x in (self.lo,self.hi) for z in (y.lo,y.hi)]
        return I(min(a)//H,ceildiv(max(a),H))
    __rmul__=__mul__
    def __truediv__(self,y):
        y=asI(y);assert y.lo*y.hi>0
        a=[Fr(x*H,z) for x in (self.lo,self.hi) for z in (y.lo,y.hi)]
        l=min(a);u=max(a)
        return I(l.numerator//l.denominator,ceildiv(u.numerator,u.denominator))
    def __rtruediv__(self,y): return asI(y)/self
    def mag(self): return max(abs(self.lo),abs(self.hi))
def asI(x): return x if isinstance(x,I) else I.exact(x)
Z=I.exact(0)
class AD:
    def __init__(self,v,d=None):
        self.v=asI(v);self.d=[Z]*4 if d is None else d
    def __add__(self,y):
        y=asAD(y);return AD(self.v+y.v,[a+b for a,b in zip(self.d,y.d)])
    __radd__=__add__
    def __neg__(self): return AD(-self.v,[-a for a in self.d])
    def __sub__(self,y): return self+-asAD(y)
    def __rsub__(self,y): return asAD(y)+-self
    def __mul__(self,y):
        y=asAD(y);return AD(self.v*y.v,[a*y.v+self.v*b for a,b in zip(self.d,y.d)])
    __rmul__=__mul__
    def __truediv__(self,y):
        y=asAD(y);v=self.v/y.v
        return AD(v,[(a-v*b)/y.v for a,b in zip(self.d,y.d)])
    def __rtruediv__(self,y): return asAD(y)/self
def asAD(x): return x if isinstance(x,AD) else AD(x)

def response(t,b):
    r=t*b;j=t+b*(1-t*t);n=j+r;T=1+b;U=1+2*b
    L=b*b+8*b+4;A=(8*b-1)*(5*b+2)
    D=5*(L*j*n-2*T*U*(j+n))
    C0=-24*D+3*A*j*n;C1=3*A*n
    C2=3*(49*b+34)*j*n-54*U*(j+n);C3=3*A*j
    B5=-4*C0+15*C1-20*C2+10*C3
    B6=-10*C0+36*C1-45*C2+20*C3
    B7=-20*C0+70*C1-84*C2+35*C3
    p=(-C0+(U+n)*B5)/D;v=(-C0+(2*U+n)*B6)/D
    a=(-C0+(T+2*j+4*r)*B7)/D;m=(-C0+(3*U+2*n)*B7)/D
    return p,v,a,m,r,D

def equations(z):
    w0=response(z[0],z[1]);w1=response(z[2],z[3])
    return [w0[0]+11,w1[0]+11,w0[1]-w1[1],w0[2]-w1[2]]
centre=[Fr(v,10**24) for v in [10504686942715102009569284,88666866113376049686925,5156906641684654809138004,167936210357531086294279]]
rho=Fr(1,10**16)
Y=[[Fr(v,10**12) for v in row] for row in [
 [74317654135440,-97353116212496,-85945879326012,20330257536319],
 [-583841728296,764803004244,675188111775,-159713860693],
 [-27730856875100,34069982371636,25456366786943,-5564302263316],
 [770570055051,-946704119133,-707367754439,154617822330]]]
box=[I(I.exact(c-rho).lo,I.exact(c+rho).hi) for c in centre]
z=[AD(v,[I.exact(int(i==j)) for j in range(4)]) for i,v in enumerate(box)]
f0=equations([I.exact(c) for c in centre]); f=equations(z)
eta=max(sum((Y[i][k]*f0[k] for k in range(4)),Z).mag() for i in range(4))
q=max(sum((I.exact(int(i==j))-sum((Y[i][k]*f[k].d[j] for k in range(4)),Z)).mag() for j in range(4)) for i in range(4))
assert eta*10**24 < H
assert q*10**6 < 5145*H
assert Fr(eta,H)+Fr(q,H)*rho<rho
bounds=[
 [(88666,88667),(931417,931418),(536056,536058),(-8215078,-8215077),(96245820,96245821),(102260466,102260467)],
 [(167936,167937),(866031,866032),(-525133,-525131),(-8215078,-8215077),(96245820,96245821),(107408353,107408354)]]
for i in range(2):
    p,v,a,m,r,D=response(box[2*i],box[2*i+1])
    for value,(lo,hi) in zip([box[2*i+1],r,D,v,a,m],bounds[i]):
        assert lo*H<value.lo*10**6 and value.hi*10**6<hi*H
if __name__=='__main__':
 print('grid_denominator',H);print('eta_upper_numerator',eta);print('q_upper_numerator',q)
 print('strict guards and all response bounds: PASS')
