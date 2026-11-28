"""各谱线对独立选取乘积与均值间的节点时，牛顿归一化后的规范投影上界。"""
from itertools import combinations_with_replacement, product
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    gram = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    source = {}
    for q, levels, direct, _ in gram["shell_weights"](root):
        nodes = [z/levels[-1] for z in levels[:-1]]+[arb(1)]
        for pivot, weights in enumerate(direct, 1):
            source[q, pivot] = list(zip(nodes, weights, strict=True))

    fields = run_path(str(root/"ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    _, lib, (central,) = fields(root, (2_345_666,), tight_neutral=True)
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    gauge = run_path(str(root/"ch39-signed-gravity-unification.py"))
    ctx.prec = 192
    pi = arb.pi()
    dot = run_path(str(root/"ch39-five-zero-gravity-gauge-transfer.py"))["dot"]
    d = (-arb(115)/218, arb(333)/218, -1)
    heavy = [(dot(d, g), z) for g, (_, _, z)
             in zip(gauge["gauge_weights"](model), central[1:], strict=True)]
    baseline = dot(d, initial)+sum(b*(arb.const_euler()+z.log()) for b, z in heavy)/(4*pi)
    multiplier = -sum(b for b, _ in heavy)/(4*pi)
    assert multiplier > 0
    pos = sum(b.max(0)*z for b, z in heavy)/(4*pi)
    neg = sum((-b).max(0)*z for b, z in heavy)/(4*pi)
    slope_bound = sum(abs(a) for _, a, _ in central)/(2*pi)

    excluded = {}
    for first, second in combinations_with_replacement(source, 2):
        nodes = [(a*b, (a+b)/2, u, v) for (a, u), (b, v)
                 in product(source[first], source[second])]
        endpoint_bounds = []
        for coherence in (-1, 1):
            absolute_first, log_upper, remainder = arb(0), arb(0), arb(0)
            for low, high, u, v in nodes:
                weight = u.real*v.real+coherence*u.imag*v.imag
                absolute_first += abs(weight)*high
                log_upper += (weight*low.log()).max(weight*high.log())
                remainder += (pos*(-weight).max(0)+neg*weight.max(0))/low
            endpoint_bounds.append((absolute_first, log_upper, remainder))
        inv_min_scale = slope_bound*max(row[0].upper() for row in endpoint_bounds)
        ceiling = baseline+max((multiplier*l+inv_min_scale*r).upper() for _, l, r in endpoint_bounds)
        assert ceiling < 0 or ceiling > 0
        if ceiling < 0:
            excluded[first, second] = ceiling
    example = excluded[(-3, 1), (-3, 2)]
    assert len(excluded) == 13 and example < arb("-4.6343")
    assert max(v.upper() for v in excluded.values()) < arb("-2.7447")
    print("给出严格负上界的配方数、最大上界、双矩构造所用配方的上界",
          len(excluded), max(v.upper() for v in excluded.values()).str(13), example.str(13))


if __name__ == "__main__":
    verify()
