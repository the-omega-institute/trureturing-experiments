"""第五壳第三枢轴仅以三个小权层承担源谱高振幅时的两组引力见证。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    number = 2_306_888_100_090
    source = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root, refine_node=number-4)[2][1]
    assert len(source) == 19 and all(v is source[0] for v in source[:17])
    m, h, k = [source[j].sqrt() for j in (0, 17, 18)]
    assert 0 < m < h < k and m*m < arb("1e-6")
    _, _, pivots, _ = next(row for row in lib("native-rh-ray-pvm")["shell_weights"](root, number) if row[0] == 6)
    c = [v.real for v in pivots[2]]
    strong = [i for i, v in enumerate(c) if abs(v) > 4]
    assert strong == [0, 9, 15, 16, 17, 18]
    moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [number])
    ctx.prec = 256
    base = moments(field, heat, m*m)
    active = (2, 3, 14)
    trace = m+h+k
    inactive = sum(v for i, v in enumerate(c) if i not in active)

    def residual(v):
        actual = [moments(field, heat, y*y) for y in v]
        f = [inactive*base[j]+sum(c[i]*r[j] for i, r in zip(active, actual, strict=True))-goal
             for j, goal in enumerate((1, arb("1e-122")))]
        return arb_mat([[f[0]], [f[1]], [sum(v)-trace]])

    def jac(v):
        d = []
        for i, y in zip(active, v, strict=True):
            ex = [(-z/y**2).exp() for _, _, z in field]
            dk = y*sum(a*t for (_, a, _), t in zip(field, ex, strict=True))/arb.pi()
            dw = -y**3*sum(g*t for (g, _, _), t in zip(field, ex, strict=True))/(2*arb.pi())
            d.append([c[i]*dk, c[i]*dw, 1])
        return arb_mat(d).transpose()

    eye = arb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    seeds = ((".2860155715730595", ".4924408711913861", ".4066318279037368"),
             (".43008477622710967", ".5279908341636202", ".22701266027745243"))
    boxes = []
    for seed in seeds:
        center = list(map(arb, seed))
        for _ in range(2):
            step = jac(center).mid().solve(residual(center).mid()).entries()
            center = [(v-a).mid() for v, a in zip(center, step, strict=True)]
        radius = arb("1e-15")
        box = [v+arb(0, radius.upper()) for v in center]
        pre = jac(center).mid().inv()
        assert not pre.det().contains(0)
        offset = max(abs(x).upper() for x in (pre*residual(center)).entries())/radius
        lip = max(sum(abs(x).upper() for x in row) for row in (eye-pre*jac(box)).tolist())
        assert offset < arb("3.51e-5") and lip < arb("1.10e-13")
        assert offset+lip < 1
        one = k-max(v.upper() for v in box)
        two = min(v.lower() for v in box)-m
        assert one > arb(".22") and two > arb(".22")
        print("第3、4、15层振幅、压缩位移比、导数误差、三维主序两项裕量",
              [x.str(16) for x in box], *[v.str(12) for v in (offset, lip, one, two)], flush=True)
        boxes.append(box)
    assert boxes[0][0] < boxes[1][0]

    # 固定另外十六层为 m，枚举除第 3、4、15 层之外的其他小权第三层。
    # 热核绝对余项 x(1-E2(z/x))、x²(1/2-E3(z/x)) 对 x 单调递增。
    pi = arb.pi()
    ka = sum(a for _, a, _ in field)/(2*pi)
    wa = -sum(g for g, _, _ in field)/(8*pi)
    xmax = k*k
    kerr = xmax*sum(abs(a)*(1-heat(2,z/xmax)) for _, a, z in field)/(2*pi)
    werr = xmax*xmax*sum(abs(g)*(arb(1)/2-heat(3,z/xmax)) for g, _, z in field)/(4*pi)
    assert ka < 0 < wa
    span = arb(0).union(k.upper())
    def disjoint(i):
        layers = (2, 3, i)
        weight = sum(abs(c[j]) for j in layers)
        rest = sum(c[j] for j in range(19) if j not in layers)
        needed = [(goal-rest*b+arb(0, (weight*e).upper()))/a for goal, b, e, a in
                  zip((1, arb("1e-122")), base, (kerr,werr), (ka,wa), strict=True)]
        def split(lo1, hi1, lo2, hi2):
            u, v = lo1.union(hi1), lo2.union(hi2)
            t = trace-u-v
            if not t.overlaps(span):
                return 1
            t = t.intersection(span)
            x = [n*n for n in (u, v, t)]
            if any(not sum(c[j]*a**p for j, a in zip(layers,x,strict=True)).overlaps(f)
                   for p, f in enumerate(needed,1)):
                return 1
            assert max(hi1-lo1,hi2-lo2)>arb("1e-6")
            if hi1-lo1 >= hi2-lo2:
                mid = (lo1+hi1)/2
                return split(lo1,mid,lo2,hi2)+split(mid,hi1,lo2,hi2)
            mid = (lo2+hi2)/2
            return split(lo1,hi1,lo2,mid)+split(lo1,hi1,mid,hi2)
        return split(arb(0), k.upper(), arb(0), k.upper())
    excluded = [(j+1, disjoint(j)) for j in range(19) if j not in (*strong,*active)]
    assert len(excluded) == 10
    print("其余第三层与排除终盒数", excluded, flush=True)
    return boxes


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
