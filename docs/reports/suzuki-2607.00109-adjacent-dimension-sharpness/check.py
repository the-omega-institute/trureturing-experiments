#!/usr/bin/env python3
"""Exact finite checks for Suzuki adjacent-dimension sharpness.

Only integers and fractions are used. These checks certify specified finite
identities and source-inequality instances, not the infinite analytic theorem.
Classical series and Machin's identity determine all interval enclosures.
"""
from fractions import Fraction as Q
from math import factorial, comb


def alternating_interval(terms):
    """Terms include one unused next term; decreasing alternating tail."""
    assert all(t >= 0 for t in terms)
    assert all(a >= b for a, b in zip(terms, terms[1:]))
    total = sum(((-1) ** k * t for k, t in enumerate(terms[:-1])), Q(0))
    other = total + (-1) ** (len(terms) - 1) * terms[-1]
    return min(total, other), max(total, other)


def atan_interval(x, n=30):
    return alternating_interval([x ** (2*k+1) / (2*k+1) for k in range(n+1)])


a, b = atan_interval(Q(1, 5)), atan_interval(Q(1, 239))
pi_lo, pi_hi = 16*a[0] - 4*b[1], 16*a[1] - 4*b[0]
assert Q(333, 106) < pi_lo < pi_hi < Q(355, 113)


def odd_df(n):
    result = 1
    for j in range(1, n+1, 2):
        result *= j
    return result


def j_integer(n, x, terms=30):
    return alternating_interval([(x/2) ** (2*k+n) / (factorial(k)*factorial(n+k))
                                 for k in range(terms+1)])


def j_half_square(n, x, terms=30):
    lo, hi = alternating_interval([x**(2*k) /
        (2**k * factorial(k) * odd_df(2*n+2*k+1)) for k in range(terms+1)])
    assert lo > 0
    return 2*x**(2*n+1)*lo**2/pi_hi, 2*x**(2*n+1)*hi**2/pi_lo


# Fidelity: source main interpolation inequality, mu=0 and nu=n+1/2.
# J_(n+1/2)(x)^2 < J_n(x)^2 + J_(n+1)(x)^2.
margins = []
for n in range(9):
    for x in [Q(1,4), Q(1,2), Q(1)]:
        jn, jnext, jhalf = j_integer(n,x), j_integer(n+1,x), j_half_square(n,x)
        assert jn[0] > 0 and jnext[0] > 0
        lower_margin = jn[0]**2 + jnext[0]**2 - jhalf[1]
        assert lower_margin > 0
        margins.append(lower_margin)
exponent = 1
while min(margins) <= Q(1, 10**exponent):
    exponent += 1
print(f'FIDELITY interpolation: 27 rational certified instances; all gaps > 1/10^{exponent}')


# Source Section 3 angular-lattice inclusion, with dimensions interpreted in R.
for d in range(2,101):
    for k in range(41):
        assert Q(k) + Q(d+2,2)-1 == Q(k+1)+Q(d,2)-1
print('FIDELITY C^(d+2) angular inclusion: d=2..100, k=0..40; 4059 identities')


# Gaussian rational arithmetic, avoiding floating-point complex numbers.
zero, one = (Q(0),Q(0)), (Q(1),Q(0))
def add(a,b): return a[0]+b[0], a[1]+b[1]
def mul(a,b): return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]
def scale(a,c): return a[0]*c, a[1]*c
def ipow(k): return [one,(Q(0),Q(1)),(Q(-1),Q(0)),(Q(0),Q(-1))][k%4]
def product(a,b):
    result={}
    for i,x in a.items():
        for j,y in b.items():
            result[i+j]=add(result.get(i+j,zero),mul(x,y))
    return result

for M in range(1,33):
    raw=product({j:one for j in range(M+1)}, {-j:one for j in range(M+1)})
    F={j:scale(v,Q(1,M+1)) for j,v in raw.items()}
    assert F[0]==one
    for j in range(-M,M+1):
        assert F[j]==(Q(M+1-abs(j),M+1),Q(0))
    assert 2*sum((Q(M+1-j,M+1) for j in range(1,M+1)),Q(0))==M
    for parity in [0,1]:
        q={j:mul(v,ipow(-j if parity else 0)) for j,v in F.items()}
        target={0:one,1:(Q(1,2),Q(0)),-1:(Q(1,2),Q(0))} if parity==0 else {
            0:one,1:(Q(0),Q(-1,2)),-1:(Q(0),Q(1,2))}
        mean=product(q,target)[0]
        assert mean==(Q(2)-Q(1,M+1),Q(0))
print('FEJER family: M=1..32; sum-of-squares expansion, mean, harmonic budget and both parity gains exact')

for d in range(2,101):
    if d%2==0:
        n=d//2-1
        n+=1-n%2
        k=n-d//2+1
        nu=Q(n)+Q(1,2)
        assert n%2==1 and k>=0
        # sin(theta - pi*(n+1/2)) = +cos(theta), since n is odd.
        assert ipow(-2*n-1)==(Q(0),Q(1))
    else:
        n=(d-1)//2
        n+=n%2
        k=n-(d-1)//2
        nu=Q(n)
        assert n%2==0 and k>=0 and ipow(-2*n)==one
    assert nu==Q(k)+Q(d+1,2)-1
for qden in range(2,101):
    delta=Q(1,qden);M=qden
    assert Q(2)-Q(1,M+1)>Q(2)-delta
print('TARGET MODES: d=2..100; nonnegative angular index and correct phase. Gain thresholds delta=1/q, q=2..100')


# Legendre polynomial coefficients, independently from the Bessel series.
def padd(a,b):
    c=a.copy()
    for i,v in b.items():c[i]=c.get(i,Q(0))+v
    return {i:v for i,v in c.items() if v}
def pscale(a,c):return {i:v*c for i,v in a.items() if v*c}
def pmul(a,b):
    c={}
    for i,v in a.items():
        for j,w in b.items():c[i+j]=c.get(i+j,Q(0))+v*w
    return {i:v for i,v in c.items() if v}
P=[{0:Q(1)},{1:Q(1)}]
for n in range(1,8):
    P.append(padd(pscale({i+1:v for i,v in P[n].items()},Q(2*n+1,n+1)),
                  pscale(P[n-1],Q(-n,n+1))))
for n in range(9):
    # Substitute t=1-u^2/2, then integrate u^(2m+1) exactly on [0,2].
    Qn={}
    for degree,c in P[n].items():
        for j in range(degree+1):
            Qn[2*j]=Qn.get(2*j,Q(0))+c*comb(degree,j)*Q(-1,2)**j
    for m in range(17):
        moment=sum((c*Q(2**(2*m+degree+2),2*m+degree+2)
                    for degree,c in Qn.items()),Q(0))
        rhs=Q((-1)**m,factorial(2*m+1))*moment
        if m<n: lhs=Q(0)
        else:
            h=m-n
            lhs=2*sum((Q((-1)**h,
                2**h*factorial(k)*factorial(h-k)*odd_df(2*n+2*k+1)*
                odd_df(2*n+2*(h-k)+1)) for k in range(h+1)),Q(0))
        assert lhs==rhs,(n,m,lhs,rhs)
print('LEGENDRE SINE IDENTITY: n=0..8; coefficients x^(2m+1), m=0..16; 153 exact matches')
print('SCOPE: finite fidelity and family identities only; no numerical quadrature, no Lean proof, no certification of the infinite sharpness theorem')
print('ALL_CHECKS_OK')
