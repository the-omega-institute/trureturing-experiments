"""One-loop 24-site threshold comparisons with fixed finite rational operations."""
from fractions import Fraction as Q


def interval(lo,hi=None):
    return Q(lo),Q(hi if hi is not None else lo)


def plus(a,b):
    return a[0]+b[0],a[1]+b[1]


def times(a,b):
    v = [x*y for x in a for y in b]
    return min(v),max(v)


def scale(x,a):
    return times(interval(x),a)


def minus(a,b):
    return plus(a,scale(-1,b))


def atan_unit(z):
    # Alternating odd powers for positive z: use a deterministic remainder bracket.
    n = 32
    s = sum((-1)**k*z**(2*k+1)/(2*k+1) for k in range(n))
    return s,s+z**(2*n+1)/(2*n+1)


def log_mantissa(m):
    assert 1<=m<=2
    z = (m-1)/(m+1)
    n = 32
    s = 2*sum(z**(2*k+1)/(2*k+1) for k in range(n))
    return s,s+2*z**(2*n+1)/((2*n+1)*(1-z*z))


def logarithm(x):
    assert x >= 1
    k = x.numerator.bit_length()-x.denominator.bit_length()
    if x<2**k:
        k -= 1
    m = x/2**k
    den = 2**96
    v = m.numerator*den//m.denominator
    lower,upper = Q(v,den),min(Q(2),Q(v+1,den))
    return plus(scale(k,log_mantissa(Q(2))), (log_mantissa(lower)[0],log_mantissa(upper)[1]))


def certify():
    pi = minus(scale(16,atan_unit(Q(1,5))),scale(4,atan_unit(Q(1,239))))
    coeff = scale(8,times(pi,pi))
    g = [Q(s) for s in (".467453",".63811","1.0583")]
    beta = [Q(41,10),-Q(19,6),-7]
    weight = [Q(1,10),Q(1,6),0]
    q0,m0 = 50017863218931,14006757711
    def t(scale_gev):
        return logarithm(Q(scale_gev,1000))
    def exposure(scale_gev,M):
        h,y = 125**2,M*M
        prod = ((h+6*y)**7*(h+4*y)**6
                *(h*h+6*h*y+4*y*y)**3*(h*h+6*h*y+6*y*y)**2)
        return scale(Q(1,2), logarithm(Q(scale_gev**46,prod)))
    def residual(i,j,T,E):
        return minus(minus(scale(1/g[i]**2-1/g[j]**2,coeff),
                           scale(beta[i]-beta[j],T)),scale(weight[i]-weight[j],E))
    def independent(scale_gev):
        return minus(scale(5/g[0]**2-3/g[1]**2-2/g[2]**2,coeff),scale(44,t(scale_gev)))
    assert independent(q0)[0]>0 and independent(q0+1)[1]<0
    T = t(q0)[0],t(q0+1)[1]
    for M,sign in [(m0,-1),(m0+1,1)]:
        E = exposure(q0,M)[0],exposure(q0+1,M)[1]
        R = residual(1,2,T,E)
        assert scale(sign,R)[0]>0
    # Eliminate radicals when checking the range in which all 23 lines are active.
    assert (3-Q(223607,100000))*m0*m0+125**2 > 12242000000**2
    assert Q(223607,100000)**2>5
    assert 6*(m0+1)**2+125**2 < 34310000000**2 < q0*q0
    for i,j in [(0,1),(1,2),(0,2)]:
        r = residual(i,j,t(q0+1),exposure(q0+1,m0))
        bound = max(abs(x) for x in r)/coeff[0]
        assert bound < Q(2,10**12)
    print("Integer GeV brackets for (Q,M)",(q0,q0+1),(m0,m0+1),
          "and all inverse-g-squared differences < 2e-12 at the stated integer inputs.")


if __name__=="__main__":
    certify()
