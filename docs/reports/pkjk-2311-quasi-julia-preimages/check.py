#!/usr/bin/env python3
"""Exact algebra and rational finite-resolution checks for arXiv:2311.13280v2.

Only integers, Fraction and exact SymPy expressions are used. The finite-resolution
basin diagnostic rounds rational coordinates to dyadics after each step; it checks
finite trajectories and does not certify membership in an infinite-time basin.
"""
from fractions import Fraction as Q
from math import isqrt
from itertools import product
import sys
import sympy as S


def zero(expr):
    assert S.cancel(S.expand(expr)) == 0, expr


def exact_checks():
    u, v, w, c, s = S.symbols('u v w c s', real=True)
    I = S.I
    rho = S.Matrix([[1+w, u-I*v], [u+I*v, 1-w]]) / 2
    H = S.Matrix([[c, s], [s, -c]])
    sq = rho.applyfunc(lambda x: x*x)
    out = H * (sq / S.trace(sq)) * H.T
    # c=cos(pi/4+epsilon/2), s=sin(pi/4+epsilon/2).
    # sin(epsilon)=s^2-c^2 and cos(epsilon)=2*s*c, c^2+s^2=1.
    coords = [out[0, 1]+out[1, 0], I*(out[0, 1]-out[1, 0]), out[0, 0]-out[1, 1]]
    expected = [(4*w*s*c+(u*u-v*v)*(s*s-c*c))/(1+w*w),
                -2*u*v/(1+w*w),
                (2*(u*u-v*v)*s*c-2*w*(s*s-c*c))/(1+w*w)]
    zero(coords[0]-expected[0]); zero(coords[2]-expected[2])
    zero(S.cancel(coords[1]-expected[1]).subs(s*s, 1-c*c))
    zero(S.trace(rho*rho)-(1+u*u+v*v+w*w)/2)
    # Source Section II: ideal pure attracting cycle and maximally mixed fixed point.
    assert F((Q(0),Q(0),Q(1)), Q(0)) == (Q(1),Q(0),Q(0))
    assert F((Q(1),Q(0),Q(0)), Q(0)) == (Q(0),Q(0),Q(1))
    assert F((Q(0),Q(0),Q(0)), Q(0)) == (Q(0),Q(0),Q(0))
    print('FIDELITY: density-matrix derivation of Eq. (8), purity, ideal cycle, C0: exact')

    t = S.symbols('t', positive=True)
    a, b = S.symbols('a b', real=True)
    eta = a+I*b
    phi = S.Matrix([[t, S.sqrt(t)*eta], [S.sqrt(t)*S.conjugate(eta), 1]])/(1+t)
    for n in range(2,7):
        power = phi.applyfunc(lambda x: x**n)
        got = power/S.trace(power)
        want = S.Matrix([[t**n, t**S.Rational(n,2)*eta**n],
                         [t**S.Rational(n,2)*S.conjugate(eta)**n, 1]])/(1+t**n)
        for x in got-want: zero(x)
        zero(got.det()-(t**n*(1-(a*a+b*b)**n))/(1+t**n)**2)
    print('IDENTITY: S_n(Phi(t,eta))=Phi(t^n,eta^n), n=2..6: symbolic, including eta=0')

    norm2 = u*u+v*v+w*w
    fp = [(2*w*c+(u*u-v*v)*s)/(1+w*w), -2*u*v/(1+w*w),
          ((u*u-v*v)*c-2*w*s)/(1+w*w)]
    defect = S.factor(1-sum(x*x for x in fp))
    defect = S.cancel(defect.subs(c*c,1-s*s))
    zero(defect-(1-norm2)*(1-w*w+u*u+v*v)/(1+w*w)**2)
    # All inverse branches for every rational interior grid point generated below.
    count = 0
    for r in [Q(-1,25),Q(0),Q(1,25)]:
        for i in range(-2,3):
            for j in range(-2,3):
                for k in range(-2,3):
                    x = (Q(i,4),Q(j,4),Q(k,4))
                    if n2(x)>=1: continue
                    y = F(x,r)
                    alt = (-x[0],-x[1],x[2])
                    assert F(alt,r)==y and n2(y)<1
                    C = trig(r)[0]*y[0]-trig(r)[1]*y[2]
                    assert C == 2*x[2]/(1+x[2]*x[2])
                    # The other quadratic root for w is outside the physical ball.
                    if x[2]: assert abs(1/x[2])>1
                    # Squaring a complex number has exactly these two roots (one at zero).
                    assert complex_sq(x)==complex_sq(alt)
                    count += 1
    print('EXACT_GRID:',count,'interior inputs; phase-independent inverse pair, physical outputs, critical axis included')


def trig(r):
    # epsilon=2 arctan(r), exact rational sine and cosine.
    return (1-r*r)/(1+r*r), 2*r/(1+r*r)


def F(x,r):
    u,v,w=x; c,s=trig(r);d=1+w*w
    return ((2*w*c+(u*u-v*v)*s)/d, -2*u*v/d,
            ((u*u-v*v)*c-2*w*s)/d)


def n2(x): return sum(t*t for t in x)
def complex_sq(x): return (x[0]*x[0]-x[1]*x[1],2*x[0]*x[1])


BITS=128
DEN=1<<BITS

def dyadic(x): return Q((x.numerator*DEN)//x.denominator,DEN)

def sqrt_approx(x):
    assert x>=0
    # Exact integer operation, returning a lower approximation with error <2^-128.
    return Q(isqrt((x.numerator << (2*BITS))//x.denominator),DEN)


def inverses_approx(y,r):
    c,s=trig(r); C=c*y[0]-s*y[2]; A=s*y[0]+c*y[2]; B=-y[1]
    assert abs(C)<1
    w=C/(1+sqrt_approx(1-C*C)); d=1+w*w
    mag=sqrt_approx(A*A+B*B)
    u=sqrt_approx(max(Q(0),(mag+A)*d/2))
    if u:
        v=B*d/(2*u)
    else:
        assert B==0
        v=sqrt_approx(max(Q(0),-A*d))
    x=tuple(map(dyadic,(u,v,w)))
    if u==0 and v==0:return [x]
    return [x,(-x[0],-x[1],x[2])]


def label(x,r):
    # Finite-time phase diagnostic only, with an even number of steps.
    for _ in range(160):x=tuple(map(dyadic,F(x,r)))
    u,v,w=x
    if abs(u)<Q(1,4) and w>Q(3,4):return 0
    if u>Q(3,4) and abs(w)<Q(1,4):return 1
    if n2(x)<Q(1,10000):return 2
    return -1


def circle(t,radius,height):
    return (radius*2*t/(1+t*t),height,radius*(1-t*t)/(1+t*t))


def finite_topology_checks():
    maps=0; boundary_cases=0; phase_pairs=0
    for n in range(1,4):
        full=(1<<n)-1; subsets=range(1<<n); topologies=[]
        for mask in range(1<<(1<<n)):
            opens={a for a in subsets if mask&(1<<a)}
            if 0 not in opens or full not in opens:continue
            if all(a|b in opens and a&b in opens for a in opens for b in opens):
                topologies.append(opens)
        assert len(topologies)==[1,4,29][n-1]
        for opens in topologies:
            def closure(a):
                closed=[full^o for o in opens if a&~(full^o)==0]
                ans=full
                for b in closed:ans&=b
                return ans
            def frontier(a):return closure(a)&closure(full^a)
            for f in product(range(n),repeat=n):
                def image(a):return sum(1<<i for i in {f[j] for j in range(n) if a&(1<<j)})
                def preimage(a):return sum(1<<j for j in range(n) if a&(1<<f[j]))
                if not all(preimage(a) in opens and image(a) in opens for a in opens):continue
                maps+=1
                for a in subsets:
                    assert preimage(frontier(a))==frontier(preimage(a))
                    boundary_cases+=1
                    for b in subsets:
                        if preimage(a)==b and preimage(b)==a:
                            assert preimage(frontier(a)&frontier(b))==frontier(a)&frontier(b)
                            assert preimage(frontier(a)|frontier(b))==frontier(a)|frontier(b)
                            phase_pairs+=1
    print('FINITE_TOPOLOGY:',maps,'continuous open self-maps on every topology of 1..3 points;',
          boundary_cases,'frontier identities;',phase_pairs,'cyclic basin pairs (intersection and union)')


def diagnostic_checks():
    # A rational negative control: endpoints in different pure phases can have a
    # third basin between them, so an arbitrary arc is not a common-boundary test.
    assert label(circle(Q(9223357,16777216),Q(99,100),Q(1,50)),Q(0))==2
    print('NEGATIVE_CONTROL: v=1/50 arc at epsilon=0 contains a mixed outcome; rejected as a two-phase bracket')
    samples=0;branches=0;max_gap=Q(0);max_residual=Q(0)
    tol=Q(1,10**6)
    for r in [Q(-1,25),Q(0),Q(1,25)]:
        for radius,height in [(Q(99,100),Q(0)),(Q(999,1000),Q(0))]:
            lo,hi=Q(0),Q(1)
            assert label(circle(lo,radius,height),r)==0
            assert label(circle(hi,radius,height),r)==1
            for _ in range(36):
                mid=(lo+hi)/2
                lab=label(circle(mid,radius,height),r)
                assert lab in [0,1],('intermediate third or unresolved basin',r,mid,lab)
                if lab==0:lo=mid
                else:hi=mid
            center=circle((lo+hi)/2,radius,height)
            left=circle(lo,radius,height);right=circle(hi,radius,height)
            ys=inverses_approx(center,r)
            ls=inverses_approx(left,r);rs=inverses_approx(right,r)
            assert len(ys)==len(ls)==len(rs)==2
            for x,xl,xr in zip(ys,ls,rs):
                assert n2(x)<1 and n2(xl)<1 and n2(xr)<1
                residual=n2(tuple(a-b for a,b in zip(F(x,r),center)))
                max_residual=max(max_residual,residual)
                assert residual<Q(1,10**60)
                gap=max(n2(tuple(a-b for a,b in zip(x,xl))),n2(tuple(a-b for a,b in zip(x,xr))))
                max_gap=max(max_gap,gap)
                assert gap<tol*tol
                # One inverse step exchanges the two phases, for each physical branch.
                assert label(xl,r)==1,(r,'left preimage',label(xl,r))
                assert label(xr,r)==0,(r,'right preimage',label(xr,r))
                branches+=1
            samples+=1
            print('BOUNDARY_DIAGNOSTIC: r='+str(r)+', radius='+str(radius)+', v='+str(height)+
                  '; 2 physical branches each bracketed by both finite-time phases within 1/10^6')
    print('DIAGNOSTICS:',samples,'boundary samples;',branches,'physical inverse branches;',
          'squared map residual < 1/10^60; squared bracketing distance < 1/10^12')
    print('LIMIT: finite-resolution phase diagnostics are not proofs of exact frontier membership or convergence')


if __name__=='__main__':
    print('Python',sys.version.split()[0],'; SymPy',S.__version__)
    exact_checks()
    finite_topology_checks()
    diagnostic_checks()
    print('ALL_CHECKS_OK')
