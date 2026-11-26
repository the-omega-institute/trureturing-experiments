"""六个时间的双高度段赔率与七节点正矩解的两组逆温度。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx


root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
gibbs = run_path(str(root/"ch39-triangle-gibbs-six-times.py"))
bands = run_path(str(root/"ch39-six-time-height-bands.py"))
rules = run_path(str(root/"ch39-seventh-time-zero-dictionary.py"))
clock = run_path(str(root/"ch39-arithmetic-height-bandwidth.py"))


def time_data():
    ctx.prec = 448
    times = clock["times"]()[:7]
    times[6] += arb("1.360900010e-8")
    d,m = bands["radial"]["compact_moments"](times)
    return times,d,m


def band_mass(lo,hi,denominator,moments):
    ends = []
    for h in (lo,hi):
        cut=1/(1+(arb(h)/rules["scale"])**2)
        nodes,weights=bands["cutoff_rule"](cut,moments,denominator)
        for upper in (False,True):
            bands["verify_indicator"](nodes,denominator,upper)
        lower=sum(v for x,v in zip(nodes,weights,strict=True) if x<cut)
        ends.append((lower,lower+weights[0]))
    return (ends[0][0]-ends[1][1]).union(ends[0][1]-ends[1][0])


def verify():
    times,d,m = time_data()
    first=band_mass(30*10**9,50*10**9,d,m)
    second=band_mass(70*10**9,100*10**9,d,m)
    ratio=second/first
    assert arb("1.5949")<ratio<arb("1.596")
    print("全谱第二段与第一段的必要算术质量比",ratio.lower().str(16),ratio.upper().str(16),flush=True)
    data=gibbs["harmonic_data"]()
    for support,title in ((rules["lower_support"],"七阶时间下端"),(rules["upper_support"],"七阶时间上端")):
        records={n:gibbs["record"](n,data,bits=112) for n in support}
        ctx.prec=448
        coordinates={n:1/(1+(h/rules["scale"])**2) for n,(h,_) in records.items()}
        _, probabilities = rules["certified_rule"](support,coordinates,times,m,d)
        w=dict(zip(support,probabilities,strict=True))
        assert all(t>0 for t in w.values())
        fitted=[]
        leverage=[]
        for lo,hi in ((30*10**9,50*10**9),(70*10**9,100*10**9)):
            sites=[n for n,(h,_) in records.items() if lo<h<hi]
            assert len(sites)==2
            i,j=sites
            h,e=records[i]
            k,f=records[j]
            log_ratio=(w[i]*h**2/(w[j]*k**2)).log()
            beta=log_ratio/(f-e)
            fitted.append(beta)
            leverage.append(abs(f-e))
            print(title, i,j,"由段内赔率反推的 β",beta.str(15),flush=True)
        assert not fitted[0].overlaps(fitted[1])
        a,b = leverage
        mismatch = abs(fitted[1]-fitted[0])*a*b/(a+b)
        assert mismatch > arb("1")
        print("任意逆温度至少需要的两点赔率倍率",mismatch.exp().str(14),flush=True)
        if support==rules["lower_support"]:
            assert arb("2.45")<fitted[0]<arb("2.47") and arb("36")<fitted[1]<arb("37")
        else:
            assert arb("1")<fitted[0]<arb("1.2") and arb("63")<fitted[1]<arb("64")


if __name__=="__main__":
    verify()
