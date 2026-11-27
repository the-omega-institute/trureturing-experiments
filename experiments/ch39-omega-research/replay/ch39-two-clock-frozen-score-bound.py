"""以源壳响应表、两钟磁谱和质量的非负性控制冻结响应表的支撑指标误差。"""
from argparse import ArgumentParser
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def scalar_mass_floor(shell, h0, h1):
    # 用 Weyl 界传播每个高度闭段，不按离散采样判断质量正性。
    scale = arb("2.435e18")/(724**2*arb("1.22089e19"))
    portal_shift = arb(".42995401872")*(arb("246.22")/arb("1.22089e19"))**2/2
    rates = [arb(0)]*19
    for i,j,length,_ in shell["geo"]["edges"]:
        if length in (724,804):
            speed = 3*abs(shell["geo"]["a" if length == 724 else "b"])
            rates[i] += speed
            rates[j] += speed
    speed = max(rates) + 12*(arb(23)/25)*arb(2)**-20/h0**3
    center = arb(15) - arb(0).union((arb(23)/25)*arb(2)**-20/h0**2/3)
    stack = [(arb(0),arb(1))]
    floors = []
    while stack:
        a,b = stack.pop()
        mid = (a+b)/2
        h = h0 + mid*(h1-h0)
        motion = speed*(b-a)*(h1-h0)/2
        levels = shell["laplacian"](3,h).eig(multiple=True)
        assert len(levels) == 19 and all(z.imag.contains(0) for z in levels)
        boxes = [z.real + arb(0,motion.upper()) for z in levels]
        mass = sorted((scale**2*v.abs_lower()**2*(v-center).abs_lower()**2).lower()
                      for v in boxes)[3]
        if mass > portal_shift:
            floors.append((mass-portal_shift).lower())
        else:
            assert b-a > arb("1e-6")
            stack.extend(((a,mid), (mid,b)))
    return min(floors), len(floors)


def common_cutoff_ceiling(root, lib, number, directions):
    # 取第一粗盒内的一点。当十九个截止相同时，所有交叉权之和为一。
    gauge = lib("all-field-signed-crossing")
    data = gauge["joint_system"](root, "1.176397e14", number=number, order=4, last=8)
    signs, masses = lib("rotation-portal-top-energy")["rotation_signs"](root, number)
    ctx.prec = 192
    mp = arb("1.22089e19")
    p0, match, probe = data["p0"], data["p1"], (arb(560000)/mp)**2
    rows, _ = lib("full-neutral-free-cutoff")["response_rows"](
        data, masses, signs, (arb("246.22")/mp)**2/2, arb(".42995401872"))
    tick = 128
    ticks = lib("fourth-third-pivot-global-kinetic-obstruction")["cutoff_ticks"]()
    assert ticks[0] < tick < 4*512
    x = p0*(arb(tick)/512).exp()
    rule = gauge["angular_rule"](13,37)
    a, b, c = data["initial"]
    initial = ((5*a/3+b-8*c/3)/12, (b-c)/3, arb(0), c/10)
    weak = lib("native-rh-ray-pvm")["weak_series"](56)(x/p0)
    # 前五行为壳无关的轻场。其余非负质量场由零质量积分包围。
    exposed = []
    for j, (kind, mass, (qq, mix, colour)) in enumerate(rows):
        terms = []
        for p in (match, probe):
            val = gauge["angular_exposure"](
                kind, mass if j < 5 else arb(0), x, p0, p, (p/p0).log(), rule)
            terms.append(val if j < 5 else arb(0).union(val))
        exposed.append((qq,mix,arb(colour),*terms))
    upper = []
    for d in directions:
        value = sum(v*q for v,q in zip(initial,d,strict=True)) + d[2]*100*weak
        for qq,mix,colour,match_kernel,probe_kernel in exposed:
            value -= ((d[0]*(qq-8*colour/3)/12 + d[1]*(mix-colour)/3)*match_kernel
                      + d[3]*colour*probe_kernel/10)/(4*arb.pi())
        upper.append((value/sum(map(abs,d))).upper())
    return max(upper)


def verify(root, *, audit=False):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    number = 2_306_888_100_097
    arc = lib("eighth-shortest-arc-pivot-support")
    support = lib("cross-shell-pivot-support")
    core = lib("fourth-third-pivot-global-kinetic-obstruction")
    tables = support["projected_tables"](core, root, number, None)
    ctx.prec = 320
    shell = lib("three-ray-chiral-boundary")
    h0,h1 = [acb.zeta_zero(n).imag for n in (number, number+1)]
    assert shell["soft"]["point"]["study"]["offset"] >= 0
    mass, leaves = scalar_mass_floor(shell, h0, h1)
    assert mass > 0
    print("全段分片数、移位后最小重标量平方质量下界", leaves, mass.str(12), flush=True)
    ceiling = common_cutoff_ceiling(root, lib, number, support["DIRECTIONS"])
    assert ceiling < arb("5.389")
    ctx.prec = 320
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)[2]
    gram_norm = sum(v*v for v in gram).sqrt()/gram[2]
    assert support["DIRECTIONS"][0] == (0,0,-1,0)
    offset, bands, scale = tables[0]
    assert offset.is_zero()
    weak_lo = arb(min(l for l,_ in bands))/scale
    weak_hi = arb(max(r for _,r in bands))/scale
    floor = (weak_lo+weak_hi-gram_norm*(weak_hi-weak_lo))/2
    assert ceiling-floor < arb("271.477")
    print("全段共同截止上界、弱差统一下界、冻结分数的误差上界",
          ceiling.str(13), floor.str(13), (ceiling-floor).str(13), flush=True)

    # 端点只读取两钟；源壳响应表也是以下所有 20³ 个格的唯一表输入。
    weights, gap = lib("eighth-ninth-wilson-height")["read_weights"](
        shell, gram, h0, shell["geo"]["phases"](h1), with_gap=True)
    # 三条电导变化各至多 |delta'(h0)| (h1-h0)；拉普拉斯对应绝对行和至多六倍。
    assert shell["soft"]["point"]["study"]["offset"] >= 0
    motion = 12*(arb(23)/25)*arb(2)**-20*(h1-h0)/h0**3
    assert gap > 2*motion
    dw = gram_norm*motion/(gap-motion)
    assert dw < arb("1e-35")
    widened = [w+arb(0,dw.upper()) for w in weights]
    major = sorted(range(19),key=lambda i: abs(weights[i]).lower(),reverse=True)[:3]
    remaining = [abs(widened[i]).upper() for i in range(19) if i not in major]
    assert min(abs(widened[i]).lower() for i in major) > max(remaining)
    total = sum(abs(w) for w in widened)
    assert total < arb("50.8")
    boxes = arc["projected_box"](tables, widened, major)
    counts = arc["SCALE"]
    lo, hi = 13283, 13287
    weak_low = boxes[0,:,:,:,0]
    weak_high = boxes[0,:,:,:,1]
    weak_floor = arb(int(weak_low.min()))/counts
    conditional_floor = int(boxes[1:,:,:,:,0].max(axis=0)[weak_low < lo].min())
    conditional_ceiling = int(boxes[1:,:,:,:,1].max(axis=0)[weak_high <= hi].min())
    # 当旧弱方向低于 lo 时，其他方向要兜住该格；上界只需要一个弱方向不高于 hi 的格。
    budget = total*arb(".004")*counts
    assert conditional_floor - budget > lo and conditional_ceiling + budget < hi
    assert tuple(map(int, boxes.max(axis=0).min(axis=(0,1,2)))) == (lo,hi)
    worst = max(arb(hi)/counts-weak_floor, ceiling-arb(lo)/counts)
    assert worst < arb("45.645")
    print("端点权误差、总变差、非弱方向两项格余量、无条件端点误差界",
          dw.str(10), total.str(12), conditional_floor, conditional_ceiling,
          worst.str(12), flush=True)

    if audit:
        target = support["projected_tables"](core, root, number+1, None)
        drift = []
        for (start, old, unit), (end, new, scale) in zip(tables, target, strict=True):
            assert unit == scale and start.overlaps(end)
            worst = max(abs(a-b) for ls,rs in zip(old,new,strict=True)
                        for a,b in zip(ls,rs,strict=True))
            drift.append(Fraction(worst,unit))
        assert drift == [0, 0, Fraction(483983,2500000000), Fraction(25707,165625000)]
        assert max(drift) < Fraction(1,250)
        print("揭晓后核对的四方向端点最大漂移", *drift, flush=True)


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("--audit", action="store_true", help="额外读取目标响应表，核对实际端点漂移")
    verify(Path(__file__).resolve().parent, audit=cli.parse_args().audit)
