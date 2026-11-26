"""固定引力截止后，最小耦合实单态的热核矩可行域和锐个数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).parent
minimal = run_path(str(ROOT / "ch39-chiral-minimal-neutral-gravity.py"))
chiral = minimal["chiral"]
E = minimal["expint"]


def required(rows,ratio):
    x = ratio**2
    vacuum = -sum(g*E(3,z/x) for g,_,z in rows)-4*arb.pi()*arb("1e-122")/x**2
    newton = 6*(2*arb.pi()/x-sum(a*E(2,z/x) for _,a,z in rows))
    assert vacuum>0 and newton>0
    return vacuum,newton


def ratio_inverse(target):
    lo,hi = arb(0),arb(64)
    assert E(2,hi)/E(3,hi)<target<2
    for _ in range(42):
        mid=(lo+hi)/2
        val=E(2,mid)/E(3,mid)
        if val>target:
            lo=mid
        elif val<target:
            hi=mid
        else:
            break
    return lo.union(hi)


def minimum_count(rows,q):
    S,T = required(rows,q)
    ratio = T/S
    assert 1<ratio<2
    # k 个零质量实标量已贡献 k 个牛顿矩。剩余质量可连续变化。
    k = int((2*S).floor().unique_fmpz())
    assert T<k   # 对以下几个尺度，最集中的质量配置的上端一定严格超过 T。
    u = ratio_inverse(ratio)
    rel = S/E(3,u)
    n = int(rel.ceil().unique_fmpz())
    assert n-1<rel<n and S<arb(n)/2
    return n,rel


def verify():
    model,_,_,_,_=chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read=chiral["source"]["spectral_supplier"](model)
    for n in (11661,11934,12207):
        rows=read(acb.zeta_zero(n).imag)
        bands={q:chiral["ordered"]([h for h,z in rows if z==q]) for q in (0,*model["charges"])}
        fields=minimal["fields"](bands)
        rejected=[]
        for scale in (arb(1),1/(8*arb.pi()).sqrt()):
            S,T=required(fields,scale)
            ratio=T/S
            assert ratio<1 if scale==1 else ratio>2
            rejected.append(ratio)
        print("零点",n,"普朗克及约化普朗克截止所需的 E₂/E₃",*[v.str(9) for v in rejected])
        for q in (arb(1)/2,1/arb(2).sqrt()):
            count, continuous=minimum_count(fields,q)
            print("Q/M_P、最少实单态数与取整前下界",q.str(6),count,continuous.str(11))


if __name__=="__main__":
    verify()
