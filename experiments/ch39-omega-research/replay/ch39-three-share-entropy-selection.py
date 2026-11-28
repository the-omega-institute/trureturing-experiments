"""在同一个已分类的三线可行集上，认证两种等权熵准则的全局唯一解。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent


def bernoulli_kl(p,q):
    if p.is_zero():
        return -(1-q).log()
    if (1-p).is_zero():
        return -q.log()
    return p*(p/q).log()+(1-p)*((1-p)/(1-q)).log()


def lower_bound(p,q):
    # Bernoulli 相对熵在 q 左侧递减，右侧递增。
    if p<q:
        return bernoulli_kl(p.upper(),q).lower()
    if p>q:
        return bernoulli_kl(p.lower(),q).lower()
    return arb(0)


def certify(point,mz,cover,reference,witness,bracket,energy_bounds):
    a,b = (arb(v) for v in bracket)
    test_energy,p,floor = point(arb(witness))
    assert all(0<v<1 for v in p) and floor>0 and mz<test_energy<arb("1.22089e19")
    limit = sum(bernoulli_kl(v,q) for v,q in zip(p,reference,strict=True)).upper()
    unit = arb(0).union(1)

    def excluded(s):
        result = point(s)
        if result is None:
            return False
        _,shares,_ = result
        if any(not unit.overlaps(v) for v in shares):
            return True
        return sum(lower_bound(v.intersection(unit),q) for v,q in zip(shares,reference,strict=True))>limit

    # 全标尺证书给出的物理标尺最小、最大端点均严格落在这里。
    outside = cover(((arb("1.00060"),a),(b,arb("1.00418"))),excluded)

    def derivatives(s):
        result = point(s,share_curvature=True)
        if result is None:
            return None
        shares,slope,curvature = result
        if not all(0<p<1 for p in shares):
            return None
        odds = [(p*(1-q)/((1-p)*q)).log() for p,q in zip(shares,reference,strict=True)]
        first = sum(d*g for d,g in zip(slope,odds,strict=True))
        second = sum(d**2/(p*(1-p))+dd*g for p,d,dd,g in zip(shares,slope,curvature,odds,strict=True))
        return first,second

    def convex(s):
        v = derivatives(s)
        return v is not None and v[1]>0

    inside = cover(((a,b),),convex)
    assert derivatives(a)[0]<0<derivatives(b)[0]
    for _ in range(45):
        mid = ((a+b)/2).mid()
        first = derivatives(mid)[0]
        if first<0:
            a = mid
        elif first>0:
            b = mid
        else:
            break
    s = a.union(b)
    energy,shares,floor = point(s)
    assert 0<min(x.lower() for x in shares) and max(x.upper() for x in shares)<1 and floor>0
    assert mz<arb(energy_bounds[0])<energy<arb(energy_bounds[1])<arb("1.22089e19")
    print("全局排除与凸性区间数、唯一标尺比、唯一能标（GeV）",
          outside,inside,s.str(12),energy.str(12),flush=True)
    return energy


def verify():
    global_energy = run_path(str(ROOT/"ch39-three-share-global-energy.py"))
    global_energy["verify"]()
    module = run_path(str(ROOT/"ch39-three-share-energy-interval.py"))
    point,mz = module["setup"]()
    repair = run_path(str(ROOT/"ch39-one-share-gravity-repair.py"))
    original = repair["source_data"]()[3]
    source = [original[i][1] for i in module["LINES"]]
    high_entropy = certify(point,mz,global_energy["cover"],[arb(1)/2]*3,"1.00365144",
                           ("1.0035","1.0038"),("8400","8405"))
    low_drift = certify(point,mz,global_energy["cover"],source,"1.00223639",
                        ("1.00215","1.00232"),("1.929e9","1.931e9"))
    assert low_drift/high_entropy>229000
    print("两个能标之比",(low_drift/high_entropy).str(12))


if __name__ == "__main__":
    verify()
