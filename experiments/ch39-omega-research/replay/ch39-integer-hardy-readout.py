"""Directed binary interval arithmetic for the finite arithmetic score on a 64-bit grid."""
from fractions import Fraction
from functools import cache
from math import comb, factorial, isqrt

UNIT = 1 << 256
GRID = 1 << 64


class Ball:
    """Closed interval with endpoints in units of 2^-256."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        self.lo, self.hi = lo, hi

    def __add__(self, other):
        return Ball(self.lo+other.lo, self.hi+other.hi)

    def __neg__(self):
        return Ball(-self.hi,-self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        values = [a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return Ball(min(values)//UNIT, -(-max(values)//UNIT))

    def scale(self, n, d=1):
        assert d>0
        a,b = sorted((n*self.lo,n*self.hi))
        return Ball(a//d,-(-b//d))

    def reciprocal(self):
        assert self.lo*self.hi > 0
        return Ball(UNIT**2//self.hi, -(-UNIT**2//self.lo))

    def square(self):
        top = max(abs(self.lo),abs(self.hi))
        low = max(0,self.lo,-self.hi)
        return Ball(low**2//UNIT,-(-top**2//UNIT))

    def sqrt(self):
        assert self.lo >= 0
        return Ball(isqrt(self.lo*UNIT),isqrt(self.hi*UNIT)+1)


def integer(v):
    return Ball(v*UNIT,v*UNIT)


ONE, ZERO = integer(1),integer(0)


def log_ratio(a,b):
    # 1 <= a/b <= 2. Bound the discarded positive atanh tail geometrically.
    assert b<=a<=2*b
    t = ONE.scale(a-b,a+b)
    t2 = t.square()
    term = t
    total = ZERO
    for j in range(88):
        total += term.scale(1,2*j+1)
        term = term*t2
    return total.scale(2) + Ball(0,term.scale(3,177).hi)


LOG_TWO = log_ratio(2,1)


@cache
def log_integer(n):
    k = n.bit_length()-1
    return LOG_TWO.scale(k) + log_ratio(n,1<<k)


def times(x,y):
    a,b=x
    c,d=y
    return (a*c-b*d, a*d+b*c)


def phase(angle):
    small = angle.scale(1,2**19)
    assert max(abs(small.lo),abs(small.hi)) <= UNIT//4
    term = ONE
    parts = [ONE,ZERO]
    for k in range(1,49):
        term = (term*small).scale(1,k)
        parts[k%2] += term.scale((-1)**(k//2))
    # Integral Taylor remainder for exp(i t) on a real argument.
    err = -(-UNIT//(4**49*factorial(49)))
    v = (parts[0]+Ball(-err,err), parts[1]+Ball(-err,err))
    for _ in range(19):
        v = times(v,v)
    return v


@cache
def arithmetic_data():
    n = 4096
    divisor = list(range(n+1))
    for p in range(2,65):
        if divisor[p]==p:
            for k in range(p*p,n+1,p):
                divisor[k]=min(divisor[k],p)
    bernoulli = [Fraction(1)]
    for m in range(1,65):
        bernoulli.append(-sum(comb(m+1,j)*x for j,x in enumerate(bernoulli))/(m+1))
    coefficients = [bernoulli[2*k]/factorial(2*k) for k in range(1,33)]
    return divisor,coefficients


def dirichlet_terms(g):
    n=4096
    divisor, _ = arithmetic_data()
    powers = [(ZERO,ZERO),(ONE,ZERO)]
    for k in range(2,n+1):
        p=divisor[k]
        if p==k:
            a,b=phase(-(g*log_integer(p)))
            factor=integer(p).sqrt().reciprocal()
            powers.append((a*factor,b*factor))
        else:
            powers.append(times(powers[p],powers[k//p]))
    return powers


def tail_coefficient(g):
    n = 4096
    _, coefficients = arithmetic_data()
    h = ONE.scale(1,2)
    inv = (g.square()+h.square()).reciprocal()
    tail = [h-(h*inv).scale(n), -(g*inv).scale(n)]
    rising = (h.scale(1,n),g.scale(1,n))
    for k,c in enumerate(coefficients,1):
        for j in (0,1):
            tail[j] += rising[j].scale(c.numerator,c.denominator)
        if k < len(coefficients):
            for shift in (2*k-1,2*k):
                rising=times(rising,((h+integer(shift)).scale(1,n),g.scale(1,n)))
    return tail


def finite_readout(g):
    powers = dirichlet_terms(g)
    n = len(powers)-1
    tail = times(powers[n],tail_coefficient(g))
    value = [sum((v[j] for v in powers[1:n]),tail[j]) for j in (0,1)]
    norm = (value[0].square()+value[1].square()).sqrt()
    # Re(2^(1/2+i g)) from the same register's n=2 term.
    return norm,powers[2][0].scale(2)


def score(m):
    g = ONE.scale(m,GRID)
    t,u = finite_readout(g)
    w = ((integer(3)-u.scale(2))*(integer(226)+u.scale(144)).sqrt()).scale(1,13)
    w = w*(g.square()+ONE.scale(1,4)).reciprocal()
    return w*(ONE+t.scale(1,256)).reciprocal()


def verify():
    # This finite shortlist is covered by the continuum-exclusion certificate.
    center = 195923358475422178854846
    scores = {m:score(m) for m in range(center-10,center+12)}
    lead = max(scores,key=lambda m:scores[m].lo)
    lower = scores[lead].lo-max(v.hi for m,v in scores.items() if m!=lead)
    assert lead==center and lower>0
    assert lower*10**31 > 175*UNIT
    assert all((v.hi-v.lo)*10**40 < UNIT for v in scores.values())
    print("Integer grid winner",lead,"absolute lead > 1.75e-29; 22 score widths < 1e-40")


if __name__ == "__main__":
    verify()
