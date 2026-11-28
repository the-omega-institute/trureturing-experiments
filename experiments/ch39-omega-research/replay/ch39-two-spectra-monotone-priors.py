"""两枚实际重谱的高能交汇所需参考概率，及共同单调谱滤波的严格上界。"""
from collections import Counter
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

ROOT = Path(__file__).resolve().parent
TARGETS = (2346136,2346144)


def load(name):
    return run_path(str(ROOT/f"ch39-{name}.py"))


def probabilities():
    model = load("three-share-energy-interval")
    cover = load("three-share-global-energy")["cover"]
    context = load("two-share-global-forward")["context"](TARGETS,20)
    *_,baseline,increments,loss = context
    assert (baseline+sum(increments[i] for i in model["LINES"]))/20+loss < 1
    bounds = [arb(".477"),arb(".261")]
    for index,heavy in enumerate(TARGETS):
        point,_ = model["setup"](heavy,batch_context=context)
        def accept(span):
            result = point(span,physical_only=True)
            if result is model["EXCLUDED"]:
                return True
            if result is None:
                return False
            energy,shares,_ = result
            if energy < arb(10)**13 or energy > arb(10)**17:
                return True
            if not all(0<p<1 for p in shares):
                return False
            _,_,_,slopes,_ = point(span,derivatives=True)
            total = sum(slopes)
            if not total>0:
                return False
            odds = sum(d*(p/(1-p)).log() for p,d in zip(shares,slopes,strict=True))/total
            q = 1/(1+(-odds).exp())
            return q>bounds[index] if index==0 else q<bounds[index]
        boxes = cover(((arb(1)/20,arb(1000)/7),),accept)
        print("目标、参考概率严格界、覆盖盒数",heavy,bounds[index].str(5),boxes,flush=True)
    return bounds[0]-bounds[1]


def verify():
    ctx.prec = 320
    exact = load("three-share-prior-obstruction")["classify"](TARGETS)
    assert exact=={2346136:[10,11],2346144:[3,4]}
    gap = probabilities()
    band = load("two-band-count-forward")
    shell = load("three-ray-chiral-boundary")
    recipes = band["verify_selection"](shell)
    spectra = load("common-spectral-band-minima")
    tables = spectra["spectra"](TARGETS)
    pairs = Counter()
    for charge,fold,intervals in recipes:
        name = "折叠谱" if fold else "原谱"
        rows = tables[charge,name,"归一化"]
        counts = tuple(band["count"](band["prefixes"](row,32),intervals) for row in rows)
        pairs[counts]+=1
    assert pairs=={(6,8):12,(6,7):8,(7,7):7,(7,8):3,(8,7):3}
    print("训练数据留下的三十三种配方在两枚新输入上的计数",sorted(pairs.items()))
    discrepancies = []
    interval_gains = []
    for key,rows in tables.items():
        prefix = low_prefix = discrepancy = gain = 0
        for _,(left,right) in spectra["level_groups"](rows):
            prefix += left-right
            discrepancy = max(discrepancy,abs(prefix))
            gain = max(gain,prefix-low_prefix)
            low_prefix = min(low_prefix,prefix)
        assert prefix==0
        discrepancies.append(discrepancy)
        interval_gains.append(gain)
        print("谱约定、绝对前缀差及正向单带差",key,discrepancy,gain)
    assert len(tables)==20 and max(discrepancies)==4 and max(interval_gains)==5
    assert gap>arb(4)/19
    print("物理所需概率差下界与单调读数差上界",gap.str(9),(arb(4)/19).str(9))


if __name__=="__main__":
    verify()
