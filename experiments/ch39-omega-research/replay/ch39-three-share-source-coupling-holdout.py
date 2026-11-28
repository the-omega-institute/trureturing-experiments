"""十二个预登记重谱迁移源端共同逆耦合，划分全部允许标尺与三份额。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent
TARGETS = tuple(range(2346100,2346123,2))


def verify():
    local = run_path(str(ROOT/"ch39-three-share-energy-interval.py"))
    domain = run_path(str(ROOT/"ch39-three-share-global-energy.py"))
    global_source = run_path(str(ROOT/"ch39-two-share-global-forward.py"))
    reference = run_path(str(ROOT/"ch39-three-share-common-coupling.py"))["source_coupling"]()
    context = global_source["context"](TARGETS,20)
    baseline,increments,loss = context[-3:]
    small = (baseline+sum(increments[i] for i in local["LINES"]))/20+loss
    assert small<arb(".769268")
    print("源端共同逆耦合、小标尺牛顿矩上界",reference.str(13),small.str(10),flush=True)
    outcomes = {}
    for number in TARGETS:
        point,mz = local["setup"](number,batch_context=context)
        roots = []
        pending = [(arb(1)/20,arb(1000)/7)]
        while pending:
            a,b = pending.pop()
            span = a.union(b)
            possible = point(span,physical_only=True,common_readout=True)
            if possible is local["EXCLUDED"] or (possible is not None and not possible[1].overlaps(reference)):
                continue
            if possible is not None:
                fit = point(span,derivatives=True,common_readout=True)
                slope = fit[3]
                if not slope.contains(0):
                    sign = 1 if slope>0 else -1
                    def residual(s):
                        return sign*(point(s,common_readout=True)[1]-reference)
                    fa,fb = residual(a),residual(b)
                    if fa>0 or fb<0:
                        continue
                    if fa<0<fb:
                        box = domain["isolate"](a,b,residual)
                        if point(box,physical_only=True) is local["EXCLUDED"]:
                            continue
                        energy,shares,floor = point(box)
                        mid = box.mid()
                        dl = point(box,derivatives=True,common_readout=True)[2]
                        energy = energy.intersection(point(mid)[0]*((box-mid)*dl).exp())
                        if any(v<0 or v>1 for v in shares) or energy<mz or energy>arb("1.22089e19"):
                            continue
                        assert all(0<v<1 for v in shares) and mz<energy<arb("1.22089e19") and floor>0
                        assert not energy.overlaps(arb(10)**13) and not energy.overlaps(arb(10)**17)
                        roots.append((box,energy))
                        continue
            assert b-a>arb("1e-10"),(number,span,possible)
            mid = ((a+b)/2).mid()
            pending.extend(((a,mid),(mid,b)))
        outcomes[number] = roots
        print(number,"全域物理解数、标尺比及交汇能标（GeV）",len(roots),
              [(s.str(10),e.str(11)) for s,e in roots],flush=True)
    expected = {2346102:("1.35178e12","1.35180e12"),
                2346108:("8.79410e11","8.79416e11"),
                2346122:("3.8682e11","3.8684e11")}
    for n,roots in outcomes.items():
        assert len(roots)==int(n in expected)
        if roots:
            lo,hi = expected[n]
            assert arb(lo)<roots[0][1]<arb(hi)<arb(10)**13
    print("物理正例", [n for n,v in outcomes.items() if v],flush=True)
    print("高能窗正例", [n for n,v in outcomes.items() if any(arb(10)**13<e<arb(10)**17 for s,e in v)],flush=True)


if __name__=="__main__":
    verify()
