"""第五、第七壳荷六解在固定谱层编号下的逐层截止距离下界。"""

from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def moment_errors(field, heat, weights, caps):
    """在各层振幅的上界处控制有质量热矩与无质量主项之差。"""
    pi = arb.pi()
    slope_k = sum(a for _, a, _ in field) / (2*pi)
    slope_w = -sum(g for g, _, _ in field) / (8*pi)
    errors = {}
    for upper in set(caps):
        x = arb(upper)**2
        error_k = x * sum(abs(a)*(1-heat(2,z/x)) for _,a,z in field) / (2*pi)
        error_w = x*x * sum(abs(g)*(arb(1)/2-heat(3,z/x)) for g,_,z in field) / (4*pi)
        errors[upper] = (error_k/abs(slope_k), error_w/abs(slope_w))
    total = [sum(abs(w)*errors[cap][j] for w, cap in zip(weights, caps, strict=True))
             for j in range(2)]
    return 1/slope_k + arb(0, total[0].upper()), arb("1e-122")/slope_w + arb(0, total[1].upper())


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    fifth, seventh = 2_306_888_100_090, 2_306_888_100_092
    _, ordered = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root, refine_node=2_306_888_100_086)[2]
    weights = []
    for node in (fifth, seventh):
        _, _, pivots, _ = next(row for row in lib("native-rh-ray-pvm")["shell_weights"](root,node)
                               if row[0] == 6)
        weights.append([v.real for v in pivots[2]])
    _, heat, fields = lib("three-index-free-atoms")["broken_fields"](root, [fifth, seventh])
    ctx.prec = 192
    cutoff_max = ".75553"
    assert ordered[-1].sqrt() < arb(cutoff_max)
    trace = sum(v.sqrt() for v in ordered)
    assert trace < arb("1.18509")
    strong = {0,9,15,16,17,18}
    assert {i for i,w in enumerate(weights[0]) if abs(w)>4} == strong
    # 第五壳六强层沿程界以及三个必需的高层分别来自已公布的独立证书。
    lower = {2:arb(".025").sqrt().lower(), 3:arb(".10").sqrt().lower(),
             14:arb(".001").sqrt().lower()}
    delta = arb(".04")
    caps5 = [".001" if i in strong else cutoff_max for i in range(19)]
    caps7 = [".041" if i in strong else cutoff_max for i in range(19)]
    second, fourth = moment_errors(fields[0],heat,weights[0],caps5)
    _, other_fourth = moment_errors(fields[1],heat,weights[1],caps7)
    u, v, w = map(arb, ("-2.25","8.45","4.1"))
    floor = u*second + v*fourth + w*other_fourth
    for i, (c, e) in enumerate(zip(*weights, strict=True)):
        lo, hi = lower.get(i,arb(0)), arb(caps5[i])
        # 取闭格覆盖第五壳振幅 y，第七壳振幅 z 在与该格相距至多 delta 的位置。
        parts = 4096
        local = None
        for j in range(parts):
            left, right = lo + (hi-lo)*j/parts, lo + (hi-lo)*(j+1)/parts
            y2 = (left**2).union(right**2)
            if e>0:
                z = min(arb(cutoff_max).upper(), (right+delta).upper())
            else:
                assert e<0
                z = max(arb(0), (left-delta).lower())
            bound = (left - c*(u*y2+v*y2**2) - e*w*z**4).lower()
            local = bound if local is None else min(local,bound)
        floor += local
    print("假定逐层差至多 0.04 时的振幅和下界、原源谱振幅和",
          floor.lower().str(15), trace.upper().str(15), flush=True)
    assert floor > arb("1.19401") > trace
    return floor


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
