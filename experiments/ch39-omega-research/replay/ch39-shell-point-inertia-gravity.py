"""Point-mass kinetic energy, ordered frequency ratios and an additive neutral threshold."""
from fractions import Fraction as F
from itertools import product
from math import isqrt
from pathlib import Path
from runpy import run_path

from flint import arb

home=Path(__file__).resolve().parent
base=run_path(str(home/"ch39-three-threshold-newton.py"))
q,A,Q,mp,baseline=[base[k] for k in ("q","A","scale","mp","baseline")]


def moments(r):
    vs=[v for v in product(range(isqrt(r)+1),repeat=3) if sum(x*x for x in v)==r]
    long=F(sum(sum(v)**2 for v in vs),3)
    trans=(len(vs)*r-long)/2
    return 2*trans/r,(trans+long)/r


def ratio(r,s,stiffness):
    ia,it=moments(r)
    ja,jt=moments(s)
    num=(stiffness/ja,stiffness/ia)
    den=(1/jt,1/it)
    return num,den


def times(a,b):
    return [a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[1]*b[1]]


def ratio_laws():
    old=ratio(17,867,F(4,9))
    new=ratio(26,29,F(8,23))
    for n,d in (old,new):
        ends=[n[i]/d[i] for i in (0,1)]
        assert all(F(4,9)<v<F(9,4) for v in ends)
        print("squared frequencies at fixed source / fixed target",*ends)
    gap=[x-y for x,y in zip(times(old[0],new[1]),times(old[1],new[0]))]
    assert all(v>0 for v in gap)
    nu,de=new
    tie=(de[0]-nu[0])/(nu[1]-de[1])
    assert tie==F(27608,29887)
    print("strict old-new gap numerator coefficients",*gap,"new isotropic when target/source mass",tie)


def main():
    ratio_laws()
    # Original two-shell alignment, all 19 target points move. Each site has inertia f^2.
    denom=q(17*867).sqrt()
    stiffness=[2*q(k)/denom for k in (396,891)]
    a,b=[k*q(1/i+1/j) for k,i,j in zip(stiffness,moments(17),moments(867))]
    groups=(1,1,2,2)
    xi=(F(1,6),F(1,3),F(7,24),F(5,24))
    assert [g*(F(1,2)-2*x) for g,x in zip(groups,xi)]==[F(v,6) for v in (1,-1,-1,1)]
    # Product of first differences of the two neutral heat translations.
    J=4*arb.pi().sqrt()/3*((a+b)**(q(3)/2)-a**(q(3)/2)-b**(q(3)/2))
    L=(a+b)*(a+b).log()-a*a.log()-b*b.log()
    assert J>0 and L>0
    def moment(k):
        return A*J*k*k.sqrt()+(-q(135)/2).union(q(19))*L*k
    def inverse(k):
        return (Q/mp)**2*(baseline+moment(k)/(12*arb.pi()))
    frac=inverse(q(1))
    print("unit axial and transverse squared masses over Q^2",a.str(14),b.str(14))
    print("predicted fraction of observed inverse Newton at unit normalization",frac.lower().str(13),frac.upper().str(13))
    target=12*arb.pi()*((mp/Q)**2-baseline)
    borders=[]
    for side in (0,1):
        lo,hi=q(0),q(10)**7
        assert moment(hi).lower()>target.upper()
        for _ in range(65):
            mid=(lo+hi)/2
            test=moment(mid).upper()<target.lower() if side==0 else moment(mid).lower()<=target.upper()
            if test:
                lo=mid
            else:
                hi=mid
        borders.append(lo if side==0 else hi)
    lo,hi=borders
    assert moment(lo).upper()<target.lower() and moment(hi).lower()>target.upper()
    # The upper envelope increases, and the lower envelope increases once positive.
    assert a>b and b>0
    print("required K/(f^2 Q^2)",lo.str(14),hi.str(14),"lightest mass /Q at least",(b*lo).sqrt().str(13))


if __name__=="__main__":
    main()
