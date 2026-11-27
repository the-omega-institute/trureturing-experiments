"""由七维源压缩的逐值谱界定量排除右外壳小转移的全部上型输运。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def support_bounds(root, delta):
    source = 2_306_888_100_086
    protocol = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    native = run_path(str(root / "ch39-free-four-site-native-positive-spectrum.py"))["native"]
    _, labels, _ = protocol["verify_nodes"](root, ((source, ".0294"),), ".1", (),
                                           ((6,7,10,13), (6,7,10,12)))
    before = native(root, source)[0]
    after = native(root, source+1)[0]
    ctx.prec = 384
    weights = labels[source,-4]
    trivial = [j for j,r in enumerate(before) if (r[0,0]+r[1,1]).real > 0]
    low, cap, high = map(arb, (".00191969", ".411496", ".947555"))
    assert len(trivial) == 7
    assert all(0 < weights[j] < cap for j in trivial)
    assert all(low < weights[j] < high for j in range(19) if j not in trivial)
    # 对源平凡块谱按其与补空间谱底的相对次序给出逐个下界。
    floor = [min(weights[j].lower(), ((1-delta)*weights[j]+delta*low).lower()) for j in trivial]
    rows = [r for r in after if (r[0,0]+r[1,1]).real > 0]
    assert len(rows) == 7
    coordinates = [((r[0,0]+r[1,1]).real, (r[0,0]-r[1,1]).real,
                    r[0,1].real, r[0,1].imag) for r in rows]
    slopes = [arb(0), arb(".000115"), arb(".0114")]
    support = [arb(1)]*len(slopes)
    for order in permutations(floor):
        tr, diff, re, im = [sum(w*z for w,z in zip(order, col, strict=True))
                            for col in zip(*coordinates, strict=True)]
        g = (tr-(diff**2+4*(re**2+im**2)).sqrt())/2
        for j,a in enumerate(slopes):
            support[j] = min(support[j], (g+a*tr).lower())
    print("平凡块三个支撑函数下界", [x.str(15) for x in support], flush=True)
    return slopes, support, (1-delta)*low, high, cap + delta*(high-cap)


def ratio_bounds(mass, slopes, support, e0, f1, upper_h):
    u,c,t = mass
    low_g = support[0]
    assert 0 < low_g < u < e0 < c < upper_h < t < f1
    # 小本征值分段，各段保留它与大本征值之间的两条仿射支撑约束。
    largest = arb(0)
    min_h = arb(1)
    bins = 1024
    for i in range(bins):
        left = low_g+(u-low_g)*i/bins
        right = low_g+(u-low_g)*(i+1)/bins
        lower_h = max(((b-(1+a)*right)/a).lower() for a,b in zip(slopes[1:],support[1:],strict=True))
        assert lower_h > c
        min_h = min(min_h, lower_h)
        ratio = (c-e0)/(e0-u)*(u-left)/(c-left)*(lower_h-u)/(lower_h-c)*(t-c)/(t-u)
        largest = max(largest, ratio.upper())
    r3 = ((t-e0)/(e0-u)*(f1-t)/(f1-u)*(u-low_g)/(t-low_g)
          * (upper_h-u)/(t-upper_h)*(t-c)/(c-u))
    assert largest < arb(".00755") and r3 < arb(".302")
    print("大本征值统一下界、两个弱轴概率比上界", min_h.str(15),
          largest.str(15), r3.str(15), flush=True)
    return ((1,arb(".00755")), (2,arb(".302")))


def cover(root, bounds):
    pencil = run_path(str(root / "ch39-triple-geometry-pencil.py"))
    lib, geom, mixing, mass, up, down, _ = pencil["setup"](root)
    ratios = ratio_bounds(mass[0], *bounds)
    normalized = [z/mass[0][-1] for z in mass[0]]
    for chart, (a,b) in enumerate(((up,down), (down,up))):
        for branch in range(3):
            stack = [(arb(-1),arb(1),0)]
            leaves = deepest = cancelled = 0
            while stack:
                lo,hi,depth = stack.pop()
                lines = lib["lines"](a,b,lo,hi,(branch,),geom)
                if branch in lines:
                    y = lines[branch]
                    prob = [abs(z)**2 for z in (mixing*y).entries()]
                    done = any(prob[0] > r*prob[j] for j,r in ratios)
                    if not done and chart == 0:
                        image = [abs(z)**2 for z in (mixing*(down*y)).entries()]
                        eigen = normalized[branch]+arb(0,max(abs(lo),abs(hi)).upper())
                        gaps = [(z-eigen)**2 for z in normalized]
                        done = any(image[0]*gaps[j] > r*image[j]*gaps[0] for j,r in ratios)
                        if done:
                            cancelled += 1
                    if done:
                        leaves += 1
                        deepest = max(deepest,depth)
                        continue
                assert depth < 26, (chart,branch,lo,hi)
                mid = (lo+hi)/2
                stack.extend(((lo,mid,depth+1),(mid,hi,depth+1)))
            print("铅笔图、支、叶区间、消参叶区间、最大深度", chart,branch,
                  leaves,cancelled,deepest,flush=True)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    cover(root, support_bounds(root, arb(".45")))
