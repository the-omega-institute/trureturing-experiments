"""固定第16、285、301线时，重谱2345712的全部可达统一能标。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent


def cover(ranges,accept):
    pending = list(ranges)
    count = 0
    while pending:
        a,b = pending.pop()
        if accept(a.union(b)):
            count += 1
            continue
        assert b-a>arb("1e-12"),(a,b)
        mid = ((a+b)/2).mid()
        pending.extend(((a,mid),(mid,b)))
    return count


def verify():
    module = run_path(str(ROOT/"ch39-three-share-energy-interval.py"))
    point,mz = module["setup"]()
    left,right = arb("1.000585"),arb("1.0042")
    outside = cover(((arb(1)/20,left),(right,arb(1000)/7)),
                    lambda s: point(s,physical_only=True) is module["EXCLUDED"])
    print("物理候选全部排除的标尺盒数",outside,flush=True)

    def regular(s):
        result = point(s,derivatives=True)
        if result is None:
            return False
        energy,shares,floor,slopes,log_slope = result
        return (0<shares[0]<1 and shares[1]<1 and 0<shares[2]<1 and floor>0
                and energy<arb("1.22089e19") and slopes[1]>0 and log_slope<0)
    inside = cover(((left,right),),regular)
    print("认证连续性、份额及能标导数的区间数",inside,flush=True)

    def isolate(lo,hi,read):
        assert read(lo)<0<read(hi)
        for _ in range(48):
            center = ((lo+hi)/2).mid()
            v = read(center)
            if v<0:
                lo = center
            elif v>0:
                hi = center
            else:
                break
        return lo.union(hi)

    scale_top = isolate(left,right,lambda s:point(s)[1][1])
    scale_ew = isolate(left,right,lambda s:mz-point(s)[0])
    assert scale_top<scale_ew
    top,shares,floor,_,log_slope = point(scale_top,derivatives=True)
    midpoint = scale_top.mid()
    center_energy = point(midpoint)[0]
    top = top.intersection(center_energy*((scale_top-midpoint)*log_slope).exp())
    assert arb("1.00060376910")<scale_top<arb("1.00060376912")
    assert arb("1.00417102485")<scale_ew<arb("1.00417102488")
    assert arb("3.095481e15")<top<arb("3.095483e15") and floor>arb(".86")
    assert arb(".58987")<shares[0]<arb(".58988") and arb(".89298")<shares[2]<arb(".893")
    print("最小、最大物理标尺比，最高能标、端点三份额、沿程最小逆耦合",
          scale_top.str(15),scale_ew.str(15),top.str(14),[p.str(11) for p in shares],floor.str(11),flush=True)


if __name__=="__main__":
    verify()
