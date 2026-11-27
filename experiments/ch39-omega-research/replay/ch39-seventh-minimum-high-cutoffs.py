"""第七壳第三枢轴的两引力条件至少需要三个高截止层。"""

from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    number = 2_306_888_100_092
    _, cutoffs = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root, refine_node=2_306_888_100_086)[2]
    source = [x.sqrt() for x in cutoffs]
    weights = next(row[2][2] for row in lib("native-rh-ray-pvm")["shell_weights"](root,number)
                   if row[0]==6)
    weights = [w.real for w in weights]
    _, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root,[number])
    ctx.prec = 192
    bounder = lib("fifth-seventh-three-high-permutations")
    bounds = bounder["thermal_bounds"](field,heat,cutoffs[-1])
    trace = sum(source)
    pair_sum = (trace-17*arb(".001")).union(trace)
    cap = source[-1].upper()
    start = (pair_sum.lower()-cap).lower()
    assert start>0

    def power(cell, exponent):
        return (cell.lower()**exponent).union(cell.upper()**exponent)

    count = 0
    for pair in combinations(range(19),2):
        c, (second, fourth) = bounder["permitted"](weights,pair,bounds,arb(0))
        stack = [(start,cap)]
        while stack:
            lo, hi = stack.pop()
            count += 1
            a = lo.union(hi)
            b = (pair_sum-a).intersection(start.union(cap))
            needed = (fourth-c[0]*power(a,4))/c[1]
            if not needed.overlaps(power(b,4)):
                continue
            b = needed.intersection(power(b,4)).sqrt().sqrt()
            if not (c[0]*power(a,2)+c[1]*power(b,2)).overlaps(second):
                continue
            assert hi-lo>arb("1e-8"),pair
            mid = (lo+hi)/2
            stack.extend(((lo,mid),(mid,hi)))
    print("第七壳第三枢轴：171 个层对的闭区间检查数",count,flush=True)
    return count


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
