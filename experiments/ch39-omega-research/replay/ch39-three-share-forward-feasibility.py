"""对十九取一检验的十二个目标放开熵限制，分类三条原线上全部标尺的高能可行性。"""
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).resolve().parent


def locate(point, excluded, lower_energy, upper_energy):
    pending = [(arb(1)/20,arb(1000)/7)]
    examined = 0
    while pending:
        lo,hi = pending.pop()
        examined += 1
        span = lo.union(hi)
        result = point(span,physical_only=True)
        if result is excluded:
            continue
        mid = ((lo+hi)/2).mid()
        if result is not None:
            energy,_,_ = result
            if energy<lower_energy or energy>upper_energy:
                continue
            witness = point(mid,physical_only=True)
            if witness is not None and witness is not excluded:
                energy,shares,floor = witness
                if lower_energy<energy<upper_energy and all(0<p<1 for p in shares) and floor>0:
                    return (mid,energy,shares,floor),examined
        assert hi-lo>arb("1e-11"),(span,result)
        pending.extend(((lo,mid),(mid,hi)))
    return None,examined


def verify():
    targets = run_path(str(ROOT/"ch39-one-of-nineteen-forward.py"))["TARGETS"]
    model = run_path(str(ROOT/"ch39-three-share-energy-interval.py"))
    context = run_path(str(ROOT/"ch39-two-share-global-forward.py"))["context"](targets,20)
    *_,baseline,increments,loss = context
    small = (baseline+sum(increments[i] for i in model["LINES"]))/20+loss
    assert small<arb(".764515")
    physical,high = [],[]
    for n in targets:
        point,mz = model["setup"](n,batch_context=context)
        witness,boxes = locate(point,model["EXCLUDED"],arb(10)**13,arb(10)**17)
        print("高能全域判定",n,witness is not None,boxes,flush=True)
        if witness is not None:
            high.append(n)
        else:
            witness,boxes = locate(point,model["EXCLUDED"],mz,arb("1.22089e19"))
            print("全物理能标判定",n,witness is not None,boxes,flush=True)
        if witness is not None:
            s,energy,shares,floor = witness
            physical.append(n)
            print("可行标尺比、能标、三份额、逆耦合下界",s.str(12),energy.str(10),
                  [p.str(8) for p in shares],floor.str(9),flush=True)
    assert physical==[2346078,2346088,2346090,2346092,2346094] and high==[2346088,2346090,2346094]
    print("全物理域中有解的编号及有高能解的编号",physical,high,flush=True)


if __name__=="__main__":
    verify()
