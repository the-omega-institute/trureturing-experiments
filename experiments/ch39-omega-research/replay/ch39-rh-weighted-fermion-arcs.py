"""Two shell arcs fix a charged row metric and a unique one-loop crossing."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

arc = run_path(str(Path(__file__).with_name("ch39-rh-internal-arc-calibration.py")))
fixed = arc["fixed"]


def exact_real(x):
    n, shift = x.man_exp()
    return Q(int(n))*Q(2)**int(shift)


def simplest(lo, hi):
    a = lo//1
    if a+1 < hi:
        return Q(a+1)
    assert a < lo < hi < a+1
    return a+1/simplest(1/(hi-a),1/(lo-a))


def ball(q):
    return arb(q.numerator)/q.denominator


def main():
    ctx.prec, ctx.cap = 1536, 3
    data = fixed["inputs"]()
    poles,rho,p,modes,k,f = data
    z, weights, _ = fixed["shell"]["quad"](21)
    def target_time(target):
        lo, hi = arb(0), arb(256)
        assert sum(w/(1+hi*v) for w,v in zip(weights,z)) < target < 1
        for _ in range(170):
            mid=(lo+hi)/2
            gap=sum(w/(1+mid*v) for w,v in zip(weights,z))-target
            if gap > 0:
                lo=mid
            elif gap < 0:
                hi=mid
            else:
                break
        return lo.union(hi)
    a,b=arb(33)/34,arb(14)/17
    first,second=target_time(a),target_time(b)
    metric=second/(3*first)
    y=(first/8).sqrt()
    print("Relative squared row norm", metric.str(16), "row multiplier", metric.sqrt().str(16))
    for bits in (16,32,64,96):
        eps=arb(2)**(-bits)
        low=target_time(b+eps)/(3*target_time(a-eps))
        high=target_time(b-eps)/(3*target_time(a+eps))
        witness=simplest(exact_real(low.upper()),exact_real(high.lower()))
        m,n=witness.numerator,witness.denominator
        assert n>1
        left_n=pow(m,-1,n)
        left_m=(m*left_n-1)//n
        left,right=Q(left_m,left_n),Q(m-left_m,n-left_n)
        assert right.numerator*left.denominator-left.numerator*right.denominator==1
        assert left.denominator+right.denominator==n
        assert ball(left)<low and low<ball(witness)<high and high<ball(right)
        # An actual finite dyadic momentum argument for this pair of integers.
        weight_ratio=3*ball(witness)
        tlo=target_time(a+eps).max(target_time(b+eps)/weight_ratio)
        thi=target_time(a-eps).min(target_time(b-eps)/weight_ratio)
        assert tlo < thi
        t=exact_real(((tlo.upper()+thi.lower())/2).mid())
        scale=2**(bits+20)
        t=Q((t*scale)//1,scale)
        assert tlo < ball(t) < thi
        for u,goal in ((ball(t),a),(ball(t)*weight_ratio,b)):
            assert abs(sum(w/(1+u*v) for v,w in zip(z,weights))-goal)<eps
        print(bits,"absolute bits: metric",m,n,"finite argument",t)
    c=[q/(fixed["CLOCK"]*y)**2 for q in poles]
    ca,cb=arb(8).log()/2,(24*metric).log()/2
    assert ca<cb
    def test(w):
        u,v=fixed["positive"](w-ca),fixed["positive"](w-cb)
        q=1000*((k-12*u-4*v)/44).exp()
        mf=q*(-w).exp()
        masses=[(mf*x.sqrt(),n) for x,n in fixed["shell"]["scalar_spectrum"](mf,c,rho,p,modes)]
        e=sum(n*fixed["positive"]((q/x).log()) for x,n in masses)
        return 22*e-402*u-46*v-15*f,q,mf,masses
    lo,hi=arb("2.48720"),arb("2.48722")
    left,right=test(lo),test(hi)
    assert left[0]<0<right[0] and right[2]*arb(8).sqrt()>1000
    active=0
    for (m,n),(v,nn) in zip(left[3],right[3]):
        assert n==nn and v>1000
        if m<right[1]:
            active+=n
        else:
            assert left[1]<v
    assert active==71 and lo>cb
    q,mf=left[1].union(right[1]),left[2].union(right[2])
    print("Weighted crossing: w",lo.union(hi).str(10),"Q",q.str(10),"M_F",mf.str(10),
          "clock unit",(mf/y).str(10),"open",active,"both fermions open",lo>cb)


if __name__=="__main__":
    main()
