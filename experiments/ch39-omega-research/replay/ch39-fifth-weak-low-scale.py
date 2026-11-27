"""第五壳的弱差归零要求三个指定层之一低于 1153 GeV。"""

from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    native = lib("native-rh-ray-pvm")
    weights = next(row[2][2] for row in native["shell_weights"](root,2_306_888_100_090)
                   if row[0]==6)
    c = [v.real for v in weights]
    _, source = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root,refine_node=2_306_888_100_086)[2]
    ctx.prec = 200
    assert [i+1 for i,v in enumerate(c) if v>0]==[1,4,10,16] and sum(c).contains(1)
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    p0 = (arb("91.1876")/arb("1.22089e19"))**2
    p1 = (arb("1.176397e14")/arb("1.22089e19"))**2
    minimum = (source[0]/p0).lower()
    threshold = (arb(1153)/arb("91.1876"))**2
    end = arb(10)**7
    weak = native["weak_series"](56)
    massive = lib("full-neutral-free-cutoff")["weak_column"]

    def cover(left,right):
        result = None
        while left<right:
            cut = min(right,(left*arb("1.003")).upper())
            cell = weak(left.union(cut))-ref
            result = cell if result is None else result.union(cell)
            left = cut
        return result

    assert 1<minimum<threshold<end
    below = cover(minimum,threshold)
    above = cover(threshold.lower(),end).union(massive(p0*end,p0,p1)-ref)
    full = below.union(above)
    # 第四层平方截止至少 0.10 时，大截止的弱核余项由单调的 1/x 界控制。
    upper = ref + (c[0]+c[9]+c[15])*above.upper() + c[3]*(massive(arb(".10"),p0,p1)-ref)
    upper += sum(w*full.lower() for w in c if w<0)
    print("三层均至少 1153 GeV 时的端点弱差上界",upper.upper().str(17),flush=True)
    assert upper < arb("-.00022")
    return upper.upper()


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
