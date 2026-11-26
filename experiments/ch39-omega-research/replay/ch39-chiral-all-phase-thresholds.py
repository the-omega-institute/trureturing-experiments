"""双钟全相位下，三代边界的十六重对与尖阈值、有限动量一圈匹配。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).parent
torus = run_path(str(ROOT / "ch39-planck-all-phase-thresholds.py"))
ordered = run_path(str(ROOT / "ch39-chiral-common-gap-crossing.py"))["ordered"]


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def sheets(model):
    bands = torus["dimensionless"](model)
    layers = {}
    for q in (0, *model["charges"]):
        layers[q] = ordered([h for h, charge, k in bands if q == charge for _ in range(k)])
        assert len(layers[q]) == (18 if q == 0 else 19)
    groups = []
    for q, weight in zip(model["charges"], model["weyl"]):
        groups.append(("fer", [2 * w for w in weight], layers[q][3:]))
    groups.append(("vec", [0, -7, -arb(21)/2], layers[0]))
    groups.append(("sca", model["higgs"], layers[3][3:]))
    groups.append(("sca", model["higgs"], layers[3][:3]))
    # 重数相同的包围只需计算一次。序号决定轻重截断，先截断再合并。
    grouped = []
    for kind, weight, masses in groups:
        hulls = {}
        for mass in masses:
            a, b = mass.lower(), mass.upper()
            assert 0 < a <= b
            key = (a.man_exp(), b.man_exp())
            if key not in hulls:
                hulls[key] = [a, b, 0]
            hulls[key][2] += 1
        grouped.append((kind, weight, list(hulls.values())))
    return grouped


def verify():
    model, _, _, initial, _ = torus["local"]["inputs"]()
    ctx.prec = 224
    groups = sheets(model)
    initial[2] = 1 / arb(".117").union(arb(".119"))
    lmin, lmax = [(arb(mu) / model["mz"]).log() for mu in ("1e13", "1e17")]
    times = lmin.union(lmax)
    # 最轻的一个重中性模对两份标量方案都存在。由 κh>M_Z 得 t<L+log h。
    assert lmax + min(b for _, b, _ in groups[5][2]).log() < 40
    pi = arb.pi()
    dirs = [(-115,333,-218), (-4,11,-7), (0,1,-1), (-1,1,0), (-1,0,1),
            (0,-1,1), (1,-1,0), (1,0,-1)]
    data = []
    for d in dirs:
        assert sum(d) == 0
        base = 2*pi*dot(d, initial) - dot(d, model["beta"])*times
        data.append((base, [dot(d, w) for _, w, _ in groups], pi*sum(map(abs,d))))
    kernels = {"sca": model["scalar"],
               "vec": lambda m, p: model["sheet"](m, p, True),
               "fer": lambda m, p: model["sheet"](m, p, False)}
    def exposures(left, right):
        # left=None 代表 t 从负无界延伸至 right。
        upper2 = (2*right).exp()
        lower2 = (2*left).exp() if left is not None else None
        ref2 = (2*(right-lmin)).exp()
        sharp, smooth = [], []
        for kind, _, masses in groups:
            total_a = total_b = arb(0)
            for a, b, k in masses:
                hi = (right-a.log()).max(0)
                lo = (left-b.log()).max(0) if left is not None else arb(0)
                total_a += k*lo.union(hi)
                above = kernels[kind](a*a,upper2)
                below = kernels[kind](b*b,lower2) if left is not None else arb(0)
                # 各核是半个归一化权重的 log(1+x) 平均，且 m≥M_Z。
                ref = (ref2/(a*a)).min(1)/8
                total_b += k*(below-ref).max(0).union(above)
            sharp.append(total_a)
            smooth.append(total_b)
        return sharp, smooth

    def lower_bounds(left, right):
        answers = []
        for exposure in exposures(left, right):
            for n in (8,7):  # 十九或十六重 Higgs。
                bounds = [((base-dot(weights[:n],exposure[:n]))/norm).lower()
                          for base, weights, norm in data]
                answers.append(max(bounds))
        return answers

    floors = [None]*4
    targets = [arb(x) for x in ("2.40","2.60","2.44","2.64")]
    covered = 0
    pending = [(arb(-8),arb(40)), (None,arb(-8))]
    while pending:
        a, b = pending.pop()
        bound = lower_bounds(a,b)
        if all(v>u for v,u in zip(bound,targets)):
            covered += 1
            floors = [v if old is None else min(old,v) for old,v in zip(floors,bound)]
        else:
            assert a is not None and b-a>arb("1e-6"), (a,b,bound)
            mid = (a+b)/2
            pending.extend(((a,mid),(mid,b)))
    print("覆盖的质量比区间数",covered)
    print("极差下界（尖/有限减法，各自十九/十六标量）",*[v.str(10) for v in floors])


if __name__ == "__main__":
    verify()
