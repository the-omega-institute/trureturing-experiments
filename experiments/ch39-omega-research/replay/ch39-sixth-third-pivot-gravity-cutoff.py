"""第六壳第三枢轴的两个引力目标至少要求一个正权层进入高截止区。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root, *, portal=None):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    number = 2_306_888_100_091
    _, ordered = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root, refine_node=2_306_888_100_086)[2]
    _, _, pivots, _ = next(row for row in lib("native-rh-ray-pvm")["shell_weights"](root, number) if row[0] == 6)
    _, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [number], portal=portal)
    ctx.prec = 192
    weights = [c.real for c in pivots[2]]
    positive = [i for i, c in enumerate(weights) if c > 0]
    negative = [i for i, c in enumerate(weights) if c < 0]
    assert [i+1 for i in positive] == [1, 10, 16] and len(negative) == 16
    pos_sum = sum(weights[i] for i in positive)
    neg_sum = -sum(weights[i] for i in negative)
    # F(x)=K(x)+250 W(x)，十九层同时达到目标时总 F 严格大于一。
    multiplier = 250
    threshold = arb(".0006")
    # 下一个二进幂的右端严格覆盖实际最大截止。
    assert ordered[-1] < 1
    edges = [arb(0), *[arb(2)**(arb(k)/16) for k in range(-480, 1)]]
    pi = arb.pi()
    previous = [(arb(0), arb(0)) for _ in field]
    floor = None
    low_ceiling = None
    for left, right in zip(edges[:-1], edges[1:], strict=True):
        current = [(right*heat(2,z/right), right**2*heat(3,z/right)) for _,_,z in field]
        # 两个正核都随截止平方非减，闭盒端点包围任意中间值。
        energy = sum(
            alpha*old_k.union(new_k)/(2*pi) - multiplier*g*old_w.union(new_w)/(4*pi)
            for (g,alpha,_), (old_k,old_w), (new_k,new_w) in zip(field, previous, current, strict=True)
        )
        floor = energy.lower() if floor is None else min(floor, energy.lower())
        if not left > threshold:
            low_ceiling = energy.upper() if low_ceiling is None else max(low_ceiling, energy.upper())
        previous = current
    bound = pos_sum*low_ceiling-neg_sum*floor
    print("单层全域下界、低截止上界与十九层线性组合上界",
          *[v.str(13) for v in (floor,low_ceiling,bound)], flush=True)
    assert bound < arb(".989853") < 1
    return tuple(i+1 for i in positive)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
