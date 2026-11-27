"""原源谱第九壳弱差为零时，第一、十六层较低能标的整 GeV 上限。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    _, squared = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root, refine_node=2_306_888_100_086)[2]
    native = lib("native-rh-ray-pvm")
    _, _, pivots, _ = next(row for row in native["shell_weights"](root, 2_306_888_100_094)
                           if row[0] == 6)
    weights = [w.real for w in pivots[2]]
    assert [i+1 for i,w in enumerate(weights) if w > 0] == [1, 10, 13, 16]
    assert all(not w.contains(0) for w in weights) and sum(weights).contains(1)

    ctx.prec = 256
    planck, mz = arb("1.22089e19"), arb("91.1876")
    p0, p1 = (mz/planck)**2, (arb("1.176397e14")/planck)**2
    reference = arb("-.074256577194663365545")
    finite = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    floor, heavy, top = (squared[i].sqrt() for i in (0, 17, 18))
    join = arb(10)**10
    tail = large(join*p0, p0, p1) - reference
    lower_correction = arb("-1.4e-11")

    def bound_from(start, upper):
        # 每个闭段调用同一带余项的有限角矩级数。最后的 1/x 包围覆盖所有更大截止。
        assert start > 1 and lower_correction < tail < upper
        cells = [(start.lower(), join)]
        count = 0
        while cells:
            left, right = cells.pop()
            response = finite(left.union(right)) - reference if right < 2*left else None
            if response is not None and lower_correction < response < upper:
                count += 1
                continue
            assert right-left > arb("1e-14")
            mid = ((left+right)/2).mid()
            cells.extend(((left, mid), (mid, right)))
        return count

    others = weights[9]+weights[12]
    important = weights[0]+weights[15]
    negative = 1-important-others
    global_upper = arb(".020483945")
    high_upper = arb(".0000272865")
    counts = (bound_from(squared[0]/p0, global_upper),
              bound_from((arb(2241)/mz)**2, high_upper))
    upper = reference + important*high_upper + others*global_upper + negative*lower_correction
    print("两条指定谱线均至少 2241 GeV 时的弱差上界、覆盖段数", upper.str(16), counts, flush=True)
    assert upper < arb("-4e-7")

    # 正权 10、13 层取源谱底；1、16 层取同一可变能标。十五条负权层均分剩余振幅迹。
    trace = 17*floor+heavy+top
    interval = arb(2240).union(arb(2241))/planck
    background = (trace-2*floor-2*interval)/15
    assert floor < interval < background < heavy < top
    # 平移 floor 后所有分量非负；总和相等且最大分量小于 top-floor。
    # 源谱只有两个非零增量，故这些条件给出全部 Schur-Horn 主序不等式。
    def example(q):
        a = arb(q)/planck
        bulk = (trace-2*floor-2*a)/15
        return (reference + important*(finite(a**2/p0)-reference)
                + others*(finite(squared[0]/p0)-reference)
                + negative*(large(bulk**2, p0, p1)-reference))

    left, right = example(2240), example(2241)
    print("同谱等能标路径在 2240 与 2241 GeV 的弱差", left.str(15), right.str(15), flush=True)
    assert left > arb(".000012") and right < arb("-4e-7")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
