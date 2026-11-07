#!/usr/bin/env python3
"""Exact rational interval checks for Haddad--Ryabogin Question 15.

Only integers and fractions.Fraction are used. Each elementary operation is
rounded outwards to a fixed rational grid. Taylor remainders are bounded by
the first omitted term; square roots are enclosed using integer square roots.
Finite checks are evidence for the proposed argument, not an infinite proof.
"""
from fractions import Fraction as Q
from math import isqrt
from itertools import product

SCALE = 10**65


def floorq(q):
    return q.numerator // q.denominator


def ceilq(q):
    return -floorq(-q)


class I:
    def __init__(self, lo, hi=None):
        lo, hi = Q(lo), Q(lo if hi is None else hi)
        assert lo <= hi
        self.lo = Q(floorq(lo * SCALE), SCALE)
        self.hi = Q(ceilq(hi * SCALE), SCALE)

    def __add__(self, other):
        o = box(other)
        return I(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -box(other)

    def __rsub__(self, other):
        return box(other) + -self

    def __mul__(self, other):
        o = box(other)
        ends = [a*b for a in (self.lo,self.hi) for b in (o.lo,o.hi)]
        return I(min(ends), max(ends))

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = box(other)
        assert not o.lo <= 0 <= o.hi, 'division through zero'
        return self * I(1/o.hi, 1/o.lo)

    def __rtruediv__(self, other):
        return box(other) / self

    def __pow__(self, n):
        assert isinstance(n,int) and n >= 0
        out = I(1)
        for _ in range(n):
            out = out * self
        return out

    def __abs__(self):
        if self.lo >= 0:
            return self
        if self.hi <= 0:
            return -self
        return I(0,max(-self.lo,self.hi))

    def short(self):
        # Outward decimal display, computed with integers, never float.
        def display(n):
            return ('-' if n < 0 else '') + str(abs(n)//10**9) + '.' + str(abs(n)%10**9).zfill(9)
        return '['+display(floorq(self.lo*10**9))+','+display(ceilq(self.hi*10**9))+']'


def box(x):
    return x if isinstance(x,I) else I(x)


def sqrt(x):
    x = box(x)
    assert x.lo >= 0
    a = isqrt(floorq(x.lo*SCALE*SCALE))
    b = isqrt(floorq(x.hi*SCALE*SCALE))
    if Q(b*b,SCALE*SCALE) < x.hi:
        b += 1
    return I(Q(a,SCALE),Q(b,SCALE))


def atan_small(x):
    x = box(x)
    assert 0 <= x.lo <= x.hi <= Q(1,2)
    n = 120
    term = x
    out = I(0)
    xx = x*x
    for j in range(n):
        out = out + (-1 if j%2 else 1)*term/(2*j+1)
        term = term*xx
    rem = term.hi/Q(2*n+1)
    return out + I(-rem,rem)


PI = 16*atan_small(Q(1,5))-4*atan_small(Q(1,239))
assert Q(314159265358979323846,10**20) < PI.lo < PI.hi < Q(314159265358979323847,10**20)


def atan(x):
    x = box(x)
    assert x.lo >= 0
    if x.hi <= Q(1,2):
        return atan_small(x)
    if x.lo > 1:
        return PI/2-atan(1/x)
    # Positive half-angle identity; new argument is < 1/2 if x <= 1.
    y = x/(1+sqrt(1+x*x))
    if y.hi <= Q(1,2):
        return 2*atan_small(y)
    return 2*atan(y)


def trig(x, sine):
    x = box(x)
    assert 0 <= x.lo <= x.hi <= 2
    n = 55
    term = x if sine else I(1)
    out = I(0)
    xx = x*x
    for j in range(n):
        out = out + (-1 if j%2 else 1)*term
        a = 2*j+2 if sine else 2*j+1
        term = term*xx/(a*(a+1))
    return out+I(-term.hi,term.hi)


def sin(x):
    return trig(x,True)


def cos(x):
    return trig(x,False)


def nearest(x):
    x = box(x)
    n = floorq(x.lo+Q(1,2))
    assert floorq(x.hi+Q(1,2)) == n
    return abs(x-n)


def delta(g,k):
    return PI*nearest((2*k-1)*g)


def ordered(g,k):
    ds = [delta(g,i) for i in range(1,k+1)]
    ds.sort(key=lambda x:x.lo)
    assert all(a.hi < b.lo for a,b in zip(ds,ds[1:]))
    return ds[0],ds[1]


def functions(g,k):
    s,t = ordered(g,k)
    R = 1/cos(PI*g)
    assert R.lo > 1 and (R*cos(s)).lo > 1
    def G(x):
        assert x.lo > 1
        return atan(sqrt(x*x-1))-atan(sqrt(R*R-1))
    # Use the literal source sine quotient, and independently its cosine form.
    ratio = sin(t-s)/(sin(t)-sin(s))
    ratio2 = cos((t-s)/2)/cos((t+s)/2)
    assert max(ratio.lo,ratio2.lo) <= min(ratio.hi,ratio2.hi)
    return G(R*ratio),-G(R*cos(s)),delta(g,k+1)


def cf_value(prefix,tail):
    x = tail
    for a in reversed(prefix):
        x = a+1/x
    return x


def cf_digits(g,n):
    x = g
    out = []
    for _ in range(n):
        a = floorq(x.lo)
        assert floorq(x.hi) == a
        out.append(a)
        x = 1/(x-a)
    return out


def convergents(digits):
    pm,p,qm,q = 1,digits[0],0,1
    out = [(p,q)]
    for a in digits[1:]:
        pm,p = p,a*p+pm
        qm,q = q,a*q+qm
        out.append((p,q))
    return out


def fidelity():
    g = (4+sqrt(2))/14
    digs = cf_digits(g,9)
    assert digs == [0,2,1,1,2,2,2,2,2]
    cs = convergents(digs)
    assert [q for _,q in cs[2:6]] == [3,5,13,31]
    for k in range(2,7):
        plus,minus,d = functions(g,k)
        assert plus.hi < d.lo and minus.hi < d.lo
    print('FIDELITY: source gamma_0=(4+sqrt(2))/14; digits and q_2..q_5=3,5,13,31; C(A_0,k) for k=2..6 verified')


def seed():
    g = (40+sqrt(2))/94
    assert cf_digits(g,12) == [0,2,3,1]+[2]*8
    assert (1/cos(PI*g)).lo > 5
    plus_bounds = [5*PI*PI/81, PI*PI/48, PI*PI/160]
    minus_bounds = [I(Q(1,2)), PI*PI/64, PI*PI/400]
    delta_bounds = [PI/5, 2*PI/25, PI/40]
    for k,pb,mb,db in zip(range(2,5),plus_bounds,minus_bounds,delta_bounds):
        plus,minus,d = functions(g,k)
        assert plus.hi < pb.lo and minus.hi < mb.lo
        assert db.hi < d.lo
        assert pb.hi < db.lo and mb.hi < db.lo
        print('SEED k='+str(k)+': Fplus='+plus.short()+', Fminus='+minus.short()+', delta_next='+d.short()+'; proposed upper/lower bounds verified')
    for k in range(5,151):
        s,t = ordered(g,k)
        plus,minus,d = functions(g,k)
        assert plus.hi < (s*t/4).lo and minus.hi < (s*s/4).lo
        assert plus.hi < d.lo and minus.hi < d.lo
    print('SEED: C(A_*,k) and stronger analytic bounds verified for k=5..150 (finite test only)')


def diophantine_small():
    cases = 0
    for ds in product(range(1,4),repeat=4):
        digits = [0,2,*ds,2,2,2]
        g = cf_value(digits[:-1],2+sqrt(2)-1)
        cs = convergents(digits)
        for j in range(3,5):
            p,q = cs[j]
            pp,qp = cs[j-1]
            odds = [(n,m) for n,m in [(q,p),(qp,pp),(q-qp,p-pp)] if n%2]
            assert len(odds) == 2 and odds[0][0] != odds[1][0]
            assert all(0 < n <= q and abs(n*g-m).hi < Q(2,q) for n,m in odds)
            D = (digits[j+1]+1)*(digits[j+2]+1)*(digits[j+3]+2)
            for k in range((q+2)//2,(cs[j+1][1]+1)//2):
                if not q <= 2*k-1 < cs[j+1][1]:
                    continue
                s,t = ordered(g,k)
                assert t.hi < (2*PI/q).lo
                assert delta(g,k+1).lo > (PI/(q*D)).hi
                cases += 1
    print('DIOPHANTINE: all '+str(cases)+' eligible (prefix,j,k) cases for a_2..a_5 in {1,2,3}, j=3,4 passed')


def tail_and_growth():
    # J=160; a_i <= i^2. Adjacent ratio of D's polynomial majorant is <2.
    def P(j):
        return ((j+1)**2+1)*((j+2)**2+1)*((j+3)**2+2)
    J = 160
    assert Q(2**(J-1)) > (4*PI.hi*P(J))**2
    # For n>=161, ((n+1)^2+c)/(n^2+c) <= (1+1/n)^2 <= (162/161)^2.
    assert Q(162,161)**12 < 2
    for n in range(1,200):
        # Integral of normalized conditional tail density on [0,1/(M+1)].
        for r in [Q(0),Q(1,3),Q(1,2),Q(1)]:
            M=n*n
            excluded=(1+r)/Q(M+1)/(1+r/Q(M+1))
            assert excluded <= Q(2,M+1)
    # For i>N, 2/(i^2+1)<2/(i(i-1)); the latter telescopes to 2/N.
    assert all(Q(2,i*i+1) < Q(2,i*(i-1)) for i in range(4,1000))
    print('TAIL/GROWTH: J=160 exponential majorant and induction ratio certified; cylinder exclusion <=2/(i^2+1); summable tail <2/N')


def main():
    fidelity()
    seed()
    diophantine_small()
    tail_and_growth()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
