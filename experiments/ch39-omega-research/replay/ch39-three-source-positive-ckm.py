"""两个荷各自任取非负谱重时，压缩三源质量锥的顶线夹角。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def source_rays(root):
    flavor = run_path(str(root / "ch39-three-gram-source-flavor.py"))
    ray = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    native = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    ctx.prec = 320
    frame, _ = flavor["family_space"](native["gram_cholesky"](root))
    adj = flavor["adj"]
    sites = [ray["soft"]["pole"], *ray["boundary_sites"]()]
    assert sites == [14,8,7,6]
    height = acb.zeta_zero(2_306_888_100_085).imag
    tables = []
    for q in (-4,2):
        matrix = ray["laplacian"](q,height)
        levels, vectors = matrix.eig(right=True,algorithm="rump")
        assert all(x.imag.contains(0) for x in levels)
        order = sorted(range(19),key=lambda k:levels[k].real.mid())
        assert all(levels[i].real < levels[j].real for i,j in zip(order,order[1:]))
        rays = []
        for k in order:
            vector = adj(frame)*acb_mat([[vectors[s,k]] for s in sites])
            norm = sum(abs(x)**2 for x in vector.entries()).sqrt()
            assert norm > 0
            rays.append(vector*(1/norm))
        tables.append(rays)
    return tables


def distance(x,y):
    # 外积给 sin(theta)，即使两线重合也不损失有效位数。
    v,w = x.entries(),y.entries()
    square = sum(abs(v[i]*w[j]-v[j]*w[i])**2 for i,j in ((0,1),(1,2),(2,0)))
    sine = arb(0).union(square.upper().sqrt()) if square.contains(0) else square.sqrt()
    cosine = abs(sum(a.conjugate()*b for a,b in zip(v,w,strict=True)))
    return arb.atan2(sine,cosine)


def sharpen(rays, index, radius, tail_ratio):
    gaps = [distance(rays[index],v) for j,v in enumerate(rays) if j != index]
    separation = min(z.lower() for z in gaps)
    assert radius < separation
    print("首层到同荷其余谱线的最小距离下界",separation.str(11))
    for _ in range(5):
        factors = [gap.sin()/(gap-radius).sin()**2 for gap in gaps]
        new = (tail_ratio*max(z.upper() for z in factors)).asin().upper()
        assert new < radius
        radius = new
    return radius


def verify(root):
    up, down = source_rays(root)
    ckm = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))
    ctx.prec = 320
    s23,s13 = [c+arb(0,ckm["sine_radius"].upper()) for c in ckm["centers"][1:3]]
    target = (s13**2 + (1-s13**2)*s23**2).sqrt().asin()
    assert arb(".0419483") < target < arb(".0420680")
    # AHS 表 2 的 1 TeV 列给各 Yukawa 的边际 1σ HPD 区间。
    published = [((".00000615",".00000014"),(".00311",".00005"),(".8616",".0043")),
                 ((".0000135",".0000002"),(".000268",".000003"),(".01401",".00009"))]
    ratios = []
    for low, middle, high in published:
        first = arb(low[0])+arb(low[1])
        second = arb(middle[0])+arb(middle[1])
        heavy = arb(high[0])-arb(high[1])
        ratios.append((first**2+second**2)/heavy**2)
    radii = [(r/(1+r)).sqrt().asin().upper() for r in ratios]
    span = sum(radii)
    print("目标重线夹角区间、初步两侧半径",target.lower().str(11),target.upper().str(11),
          *(r.str(11) for r in radii))
    small, large, pending = [], [], []
    for i, u in enumerate(up):
        small_part = []
        for j, d in enumerate(down):
            angle = distance(u,d)
            assert angle.is_finite()
            if angle+span < target:
                small.append(angle)
                small_part.append(j+1)
            elif angle-span > target:
                large.append(angle)
            else:
                pending.append((i,j,angle))
        print("上型层，初筛角度过小的下型层",i+1,small_part)
    assert (len(small),len(large),len(pending)) == (73,287,1)
    [(i,j,angle)] = pending
    assert (i,j) == (0,0) and arb(".03726589") < angle < arb(".037265892")
    assert max(z.upper() for z in small) < arb(".017748728")
    assert min(z.lower() for z in large) > arb(".136020310")
    up_radius = sharpen(up,i,radii[0],ratios[0])
    down_radius = sharpen(down,j,radii[1],ratios[1])
    assert up_radius < arb(".00006664") and down_radius < arb(".00389692")
    small_max = max((max(z.upper() for z in small)+span).upper(),
                    (angle+up_radius+down_radius).upper())
    large_min = (min(z.lower() for z in large)-span).lower()
    assert small_max < arb(".04122946") < target < arb(".1128435") < large_min
    assert small_max.cos() > arb(".99915018") and large_min.cos() < arb(".9936401")
    # 同表 CKM 的输入为弧度。余弦在两个正角区间上均递减。
    table_low = arb(".04317").cos()*arb(".00385").cos()
    table_high = arb(".04233").cos()*arb(".00369").cos()
    assert large_min.cos() < arb(".99906091") < table_low < table_high
    assert table_high < arb(".99909742") < small_max.cos()
    print("表 2 同尺度 CKM 顶底模的两端",table_low.str(13),table_high.str(13))
    print("最终两侧半径和两支角度的界",up_radius.str(11),down_radius.str(11),
          small_max.str(11),large_min.str(11))
    print("全域可行角度两支对应的模界，参考模上、下端",
          small_max.cos().str(11),large_min.cos().str(11),
          target.lower().cos().str(11),target.upper().cos().str(11))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
