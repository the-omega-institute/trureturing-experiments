"""预登记的十二枚重谱，固定十九取一伯努利参考概率的内点预报。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
TARGETS = tuple(range(2346076,2346099,2))


def derivative_bound(p, d, z, direction):
    """给 direction*L' 的上界。返回 None 表示区间的一端无界。"""
    total = arb(0)
    for v,dp in zip(p,d,strict=True):
        slope = direction*dp
        if v>0 and v<1:
            term = slope*((v/(1-v)).log()-z)
        elif slope>0 and v.upper()<1 and v.upper()>0:
            h = v.upper()
            term = slope*((h/(1-h)).log()-z)
        elif slope<0 and v.lower()>0 and v.lower()<1:
            h = v.lower()
            term = slope*((h/(1-h)).log()-z)
        else:
            return None
        total += term.upper()
    return total


def search(point,z, *, exclude=None):
    def first(s):
        _,p,_,dp,_ = point(s,derivatives=True)
        assert all(0<x<1 for x in p)
        return sum(d*((v/(1-v)).log()-z) for v,d in zip(p,dp,strict=True))

    queue = [(arb(".95"),arb("1.05"))]
    roots = []
    boxes = 0
    while queue:
        lo,hi = queue.pop()
        s = lo.union(hi)
        boxes += 1
        if exclude is not None and exclude(s):
            continue
        result = point(s,share_curvature=True)
        if result is not None:
            p,dp,ddp = result
            if any(v<0 or v>1 for v in p):
                continue
            upper = derivative_bound(p,dp,z,1)
            negative_upper = derivative_bound(p,dp,z,-1)
            if (upper is not None and upper<0) or (negative_upper is not None and negative_upper<0):
                continue
            if all(0<v<1 for v in p):
                second = sum(d*d/(v*(1-v))+dd*((v/(1-v)).log()-z)
                             for v,d,dd in zip(p,dp,ddp,strict=True))
                if second>0:
                    left,right = first(lo),first(hi)
                    if left>0 or right<0:
                        continue
                    if left<0<right:
                        for _ in range(45):
                            mid = ((lo+hi)/2).mid()
                            val = first(mid)
                            if val<0:
                                lo = mid
                            elif val>0:
                                hi = mid
                            else:
                                break
                        roots.append(lo.union(hi))
                        continue
        assert hi-lo>arb("1e-12"),(s,result)
        mid = ((lo+hi)/2).mid()
        queue.extend(((lo,mid),(mid,hi)))
    return roots,boxes


def verify():
    load = lambda name: run_path(str(ROOT/f"ch39-{name}.py"))
    setup = load("three-share-energy-interval")["setup"]
    shared = load("two-share-global-forward")["context"](TARGETS,20)
    model,_,_,_,_ = load("chiral-minimal-neutral-gravity")["chiral"]["source"]["mean"]["inputs"]()
    beta = [v/(2*arb.pi()) for v in model["beta"]]
    z = -arb(18).log()
    physical,high,stationary = [],[],[]
    for heavy in TARGETS:
        point,mz = setup(heavy,batch_context=shared)
        roots,boxes = search(point,z)
        assert len(roots)<=1
        if roots:
            stationary.append(heavy)
        print("驻点数及标尺盒数",heavy,len(roots),boxes,flush=True)
        for root in roots:
            energy,common = point(root,common_readout=True)
            _,shares,_ = point(root)
            ell = (energy/mz).log()
            couplings = [common,*[common+b*ell for b in beta]]
            assert all(not v.contains(0) for v in couplings) and all(0<p<1 for p in shares)
            predicates = [ell,(arb("1.22089e19")/energy).log(),*couplings]
            assert all(not v.contains(0) for v in predicates)
            good = all(v>0 for v in predicates)
            window = good and arb(10)**13<energy<arb(10)**17
            if good:
                assert not (energy-arb(10)**13).contains(0) and not (energy-arb(10)**17).contains(0)
                physical.append(heavy)
            if window:
                high.append(heavy)
            labels = ("低于电弱","高于普朗克","共同逆耦合","第一支逆耦合","第二支逆耦合","第三支逆耦合")
            failed = [label for label,v in zip(labels,predicates,strict=True) if v<0]
            print("目标、标尺比、能标、物理与高能判定、三概率、越界条件",heavy,root.str(12),energy.str(10),
                  good,window,[p.str(8) for p in shares],failed,flush=True)
    assert stationary==[2346076,2346078,2346088,2346090,2346092,2346094]
    assert physical==[2346078,2346092] and not high
    print("物理正例与书中高能窗",physical,high,flush=True)


if __name__ == "__main__":
    verify()
