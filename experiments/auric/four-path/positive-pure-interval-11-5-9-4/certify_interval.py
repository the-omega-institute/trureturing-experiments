"""Exact closed-interval certificates for the Section47 interval consumer.

Only Python standard-library integer/rational arithmetic is used for proof.
All square roots are rounded outward to multiples of 2**(-96).
Sources: 5fd6a884f9173f778cfa040217d7d51f748ac9c6, FourPath
Continuation I (37.11)--(37.22), Continuation II (46.11)--(46.13),
(47.1), and (47.9)--(47.13). Output is mathematical certificate data.
"""
from fractions import Fraction as Q
from math import isqrt
import hashlib
import json
from pathlib import Path

BITS = 96
SCALE = 1 << BITS
OUT = Path(__file__).resolve().with_name('certificate.json')

class IV:
    __slots__ = ('lo', 'hi')
    def __init__(self, lo, hi=None):
        self.lo = Q(lo)
        self.hi = self.lo if hi is None else Q(hi)
        assert self.lo <= self.hi
    @staticmethod
    def of(v):
        return v if isinstance(v, IV) else IV(v)
    def __add__(self, v):
        v = IV.of(v)
        return IV(self.lo + v.lo, self.hi + v.hi)
    __radd__ = __add__
    def __neg__(self):
        return IV(-self.hi, -self.lo)
    def __sub__(self, v):
        return self + (-IV.of(v))
    def __rsub__(self, v):
        return IV.of(v) - self
    def __mul__(self, v):
        v = IV.of(v)
        a = [self.lo*v.lo, self.lo*v.hi, self.hi*v.lo, self.hi*v.hi]
        return IV(min(a), max(a))
    __rmul__ = __mul__
    def inv(self):
        assert self.lo > 0 or self.hi < 0, 'zero in divisor'
        return IV(1/self.hi, 1/self.lo)
    def __truediv__(self, v):
        return self*IV.of(v).inv()
    def __rtruediv__(self, v):
        return IV.of(v)*self.inv()
    def sq(self):
        a = self.lo*self.lo
        b = self.hi*self.hi
        return IV(0 if self.lo <= 0 <= self.hi else min(a,b), max(a,b))
    def sqrt(self):
        assert self.lo >= 0, 'negative radicand'
        n = isqrt((self.lo.numerator*SCALE*SCALE)//self.lo.denominator)
        m = isqrt((self.hi.numerator*SCALE*SCALE)//self.hi.denominator)
        lo = Q(n,SCALE)
        hi = Q(m,SCALE) if Q(m*m,SCALE*SCALE) == self.hi else Q(m+1,SCALE)
        assert lo*lo <= self.lo and self.hi <= hi*hi
        return IV(lo,hi)
    def serial(self):
        return [str(self.lo),str(self.hi)]

h = IV(Q(1,2)).sqrt()
c = ((1+h)/2).sqrt()
s = ((1-h)/2).sqrt()
k = c-h
s2 = (1-h)/2

def omega(xi,q):
    xi,q = IV.of(xi),IV.of(q)
    p = 1-q
    x = 1-xi/4+q/2
    y = (1-x.sq()).sqrt()
    R0 = h*x+c
    mu = (h*k*p.sqrt()/R0+h*q.sqrt()/y)/2
    nu = (k*p.sqrt()+(IV(Q(1,2))*x-s2)*q.sqrt()/(h*y))/2
    tc2 = s2*q*(1+x)/x.sq()
    return nu-tc2/(4*mu)

def derivative(xi,q=Q(9,20)):
    xi,q = IV.of(xi),IV.of(q)
    p = 1-q
    x = 1-xi/4+q/2
    y = (1-x.sq()).sqrt()
    R0 = h*x+c
    gam = p.sqrt()*R0+q.sqrt()*h*y
    bracket = h*p.sqrt()-R0/p.sqrt()+h*y/q.sqrt()-h*x*q.sqrt()/y
    return -(4*h+8*s)-4*gam*bracket

def floorQ(v,n=10**9):
    return Q((v.numerator*n)//v.denominator,n)

def ceilQ(v,n=10**9):
    return -floorQ(-v,n)

def certify(name,fn,lo,hi,sense,target):
    leaves = []
    def visit(a,b,depth):
        value = fn(IV(a,b))
        passes = value.lo >= target if sense == '>=' else value.hi <= target
        if passes:
            leaves.append({'a':str(a),'b':str(b),'value':value.serial()})
            return
        assert depth < 24, (name,a,b,value.serial())
        m=(a+b)/2
        visit(a,m,depth+1)
        visit(m,b,depth+1)
    visit(lo,hi,0)
    # The accepted closed intervals partition the original closed domain:
    assert Q(leaves[0]['a'])==lo and Q(leaves[-1]['b'])==hi
    assert all(Q(t['a'])<Q(t['b']) for t in leaves)
    assert all(Q(a['b'])==Q(b['a']) for a,b in zip(leaves,leaves[1:]))
    assert sum((Q(t['b'])-Q(t['a']) for t in leaves),Q(0))==hi-lo
    lower=min(Q(t['value'][0]) for t in leaves)
    upper=max(Q(t['value'][1]) for t in leaves)
    return {'name':name,'domain':[str(lo),str(hi)],'sense':sense,
            'target':str(target),'pieces':len(leaves),
            'aggregate_enclosure':[str(floorQ(lower)),str(ceilQ(upper))],
            'leaves':leaves}

def main():
    # All these are mathematical checks, not sampled sign checks.
    certificates = [
        certify('C1',lambda q:omega(Q(9,4),q),Q(0),Q(1,10),'>=',Q(1,10000)),
        certify('C2',lambda q:omega(Q(9,4),q),Q(3,8),Q(1),'>=',Q(1,10000)),
        certify('C3',lambda q:omega(Q(9,4),q),Q(1,10),Q(3,8),'>=',-Q(21,1000)),
        certify('C4',derivative,Q(11,5),Q(9,4),'<=',-Q(1,10)),
    ]
    point = omega(Q(9,4),Q(1,4))
    assert point.hi < -Q(19,1000)
    # Universal constants used in the analytic continuous-domain deduction.
    assert (8*h*c-1).lo > 4
    assert s.lo > Q(19,50)
    assert k.lo > Q(21,100)
    assert (h*Q(51,80)+c).hi < Q(7,5)
    assert (h+c).hi < Q(5,3)
    assert IV(30).sqrt().hi < Q(11,2)
    b_floor=Q(38,275)+Q(6615,2240000)
    assert b_floor > Q(7,50)
    gain = Q(21,1000)**2/(4*Q(7,50))
    assert gain == Q(63,80000)
    separation = Q(9,800)-Q(40,3)*gain-4*gain**2
    assert separation == Q(1196031,1600000000) > Q(7,10000)
    # Whole pure I-domain: q in [0,1], x in [7/16,19/20].
    domain_x=1-IV(Q(11,5),Q(9,4))/4+IV(0,1)/2
    assert domain_x.lo == Q(7,16) and domain_x.hi == Q(19,20)
    sign_gap=18-8*(c+s)-2*Q(9,4)
    assert sign_gap.lo > 3
    result={
        'source_pin':'5fd6a884f9173f778cfa040217d7d51f748ac9c6',
        'arithmetic':'Fraction; isqrt outward rounding; 96 binary fractional bits',
        'rule':'Natural interval extension, recursive closed midpoint bisection; '
               'every leaf satisfies the named rational inequality; exact contiguous coverage checks.',
        'certificates':certificates,
        'omega_9_4_at_1_4':point.serial(),
        'constants':{'h':h.serial(),'c':c.serial(),'s':s.serial(),
                     'curvature':(8*h*c-1).serial(),
                     'B_strict_lower':str(b_floor),
                     'J_upper':str(gain),'strict_gap_lower':str(separation),
                     'gap_excess_over_7_10000':str(separation-Q(7,10000)),
                     'x_domain':domain_x.serial(),
                     'inactive_sign_gap_lower':str(floorQ(sign_gap.lo))},
    }
    OUT.write_text(json.dumps(result,indent=2)+'\n')
    data=OUT.read_bytes()
    print(json.dumps({'status':'PASS','certificates':[
        {a:b for a,b in t.items() if a!='leaves'} for t in certificates],
        'point_omega_enclosure':[str(floorQ(point.lo)),str(ceilQ(point.hi))],
        'strict_gap_lower':str(separation),
        'excess_over_7_10000':str(separation-Q(7,10000)),
        'coverage_manifest':str(OUT),'manifest_sha256':hashlib.sha256(data).hexdigest()},indent=2))

if __name__=='__main__':
    main()
