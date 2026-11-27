"""第五壳第三枢轴的五式、源谱主序与中性规范路径。"""

from inspect import getclosurevars
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def solve(root, centers, axes, *, number=2_306_888_100_090, pivot=3, tighten=False):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    _, ordered = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root, refine_node=2_306_888_100_086
    )[2]
    source = [x.sqrt() for x in ordered]
    native = lib("native-rh-ray-pvm")
    _, _, pivots, _ = next(
        a for a in native["shell_weights"](root, number) if a[0] == 6
    )
    assert 1 <= pivot <= len(pivots)
    weights = [w.real for w in pivots[pivot - 1]]
    assert sum(weights).contains(1)
    gauge = lib("all-field-signed-crossing")
    neutral = lib("full-neutral-free-cutoff")
    data = gauge["joint_system"](root, "1.176397e14", number=number, order=4, last=8)
    signs, masses = lib("rotation-portal-top-energy")["rotation_signs"](root, number)
    moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](
        root, [number]
    )
    ctx.prec = 192
    rho = (arb("246.22") / arb("1.22089e19")) ** 2 / 2
    species, _ = neutral["response_rows"](data, masses, signs, rho, arb(".42995401872"))
    p0, p1 = data["p0"], data["p1"]
    pi = arb.pi()
    log = (p1 / p0).log()
    a1, a2, a3 = data["initial"]
    terms = [
        (k, z, ((q - arb(8) * c / 3) / (4 * pi), (m - c) / (4 * pi)))
        for k, z, (q, m, c) in species
    ]
    means = lib("fourth-third-pivot-cyclic-five")["angular_means"]
    rule = gauge["angular_rule"](17, 46)
    weak = native["weak_series"](56)
    poly = getclosurevars(weak).nonlocals["coefficients"]

    def weak_at(x):
        return weak(x / p0) if x < 10**7 * p0 else neutral["weak_column"](x, p0, p1)

    def columns(x):
        cols = [arb(0), arb(0)]
        for k, z, coeff in terms:
            v = gauge["angular_exposure"](k, z, x, p0, p1, log, rule)
            for j in range(2):
                cols[j] -= coeff[j] * v
        return [*cols, weak_at(x), *moments(field, heat, x)]

    def jac(sites):
        cols = []
        for i in axes:
            x = sites[i]
            ang = [a - b for a, b in zip(means(p0 / x), means(p1 / x), strict=True)]
            d = [arb(0), arb(0)]
            for k, z, b in terms:
                v = (-z / x).exp() * ang[k]
                for j in range(2):
                    d[j] -= b[j] * v
            dw = (
                -sum(n * a * (p0 / x) ** n for n, a in enumerate(poly, 1))
                if x < 10**7 * p0
                else arb(0)
            )
            e = [(-z / x).exp() for _, _, z in field]
            dk = (
                x * sum(a * v for (_, a, _), v in zip(field, e, strict=True)) / (2 * pi)
            )
            dv = (
                -x
                * x
                * sum(g * v for (g, _, _), v in zip(field, e, strict=True))
                / (4 * pi)
            )
            cols.append([weights[i] * a for a in (*d, dw, dk, dv)] + [x.sqrt() / 2])
        return arb_mat(cols).transpose()

    sites = list(map(arb, centers))
    readings = [columns(x) for x in sites]

    def residual():
        f = [
            a + sum(c * row[j] for c, row in zip(weights, readings, strict=True))
            for j, a in enumerate(
                (5 * a1 / 3 + a2 - 8 * a3 / 3, a2 - a3, 0, -1, -arb("1e-122"))
            )
        ]
        return arb_mat([[a] for a in [*f, sum(x.sqrt() for x in sites) - sum(source)]])

    # 数值中心的修正在同一组六个坐标上完成。最终证书重新计算残差。
    for _ in range(1):
        shift = jac(sites).mid().solve(residual().mid()).entries()
        for i, d in zip(axes, shift, strict=True):
            sites[i] = (sites[i] * (-d).exp()).mid()
            readings[i] = columns(sites[i])
    radius = arb("2e-7")
    box = [
        x * arb(0, radius.upper()).exp() if i in axes else x
        for i, x in enumerate(sites)
    ]
    pre = jac(sites).mid().inv().mid()
    assert not pre.det().contains(0)
    f = residual()
    ref = arb("-.074256577194663365545")
    # 对真实弱差和光滑替身间的差，在盒内与中心各留一次余量。
    err = arb(0)
    for i in axes:
        x = box[i]
        if x < 10**7 * p0:
            err += abs(weights[i]) * (arb("1e-19") + 2 * x / p1)
        else:
            err += abs(weights[i]) * abs(weak_at(x) - ref)
    f[2, 0] += arb(0, (2 * err).upper())
    eye = arb_mat([[int(i == j) for j in range(6)] for i in range(6)])
    offset = (pre * f).entries()
    mismatch = (eye - pre * jac(box)).tolist()
    ratio = max(
        (abs(v) / radius + sum(map(abs, row))).upper()
        for v, row in zip(offset, mismatch, strict=True)
    )
    print("六维自映比", ratio.str(13), flush=True)
    assert ratio < 1
    if tighten:
        # 真根偏移逐坐标满足预条件后的非负均值不等式。封闭盒中原先已覆盖真实弱差与平滑替身的差。
        widths=[radius]*len(axes)
        for _ in range(10):
            mismatch=(eye-pre*jac(box)).tolist()
            bounds=[abs(v).upper()+sum(abs(a).upper()*w for a,w in zip(row,widths,strict=True))
                    for v,row in zip(offset,mismatch,strict=True)]
            widths=[min(old,new.upper()) for old,new in zip(widths,bounds,strict=True)]
            for i,w in zip(axes,widths,strict=True):
                box[i]=sites[i]*arb(0,w.upper()).exp()
        print("原五式各调整坐标的收紧半宽",*[v.str(9) for v in widths],flush=True)
    gaps = lib("schur-horn-gravity-holdout")["ky_fan_gaps"](
        source, [x.sqrt() for x in box]
    )
    print("源谱主序最小余量", min(gaps).str(14), flush=True)
    assert min(gaps) > 0
    w, z, t = [
        k * rho
        for k in (
            arb(".63811") ** 2 / 2,
            (arb(".63811") ** 2 + 3 * arb(".467453") ** 2 / 5) / 2,
            arb(".8616") ** 2,
        )
    ]
    coeff = (
        (2 * w + z) / (24 * pi) + (w + z) / (48 * pi) + t / (16 * pi) + t / (24 * pi)
    )
    # 按逆截止对有符号的光滑余项分部求和。相等的固定位置产生零长度增量。
    last = arb(0)
    total = arb(1)
    variation = arb(0)
    order = sorted(range(19), key=lambda i: sites[i].mid(), reverse=True)
    previous = None
    for i in order:
        inv = 1 / box[i]
        if not inv > last:
            assert (
                previous not in axes
                and i not in axes
                and centers[i] == centers[previous]
            )
        variation += abs(inv - last) * abs(total)
        total -= weights[i]
        last = inv
        previous = i
    weak_bound = (
        lib("full-neutral-running")["weak_bound"](
            data, [p0], [arb(0)], lib("neutral-axial-top"), lib("neutral-bosonic")
        )
        + coeff * variation
    )
    print("全程弱差绝对界", weak_bound.upper().str(13), flush=True)
    return data, species, weights, sites, box, weak_bound, lib


def verify(root):
    data, species, weights, sites, box, weak, lib = solve(
        root,
        (
            "4.59881149198487819e-33",
            "1.98703768478281753e-34",
            "0.081804926717849849",
            "0.242497971196136003",
            "2.01802914679143173e-34",
            "1.64094530819905405e-15",
            "7.42156679413875031e-19",
            "2.155180747019892e-34",
            "1.65647727903706752e-34",
            "7.12319811999282302e-35",
            "2.69670158759411499e-34",
            "7.76446987057365487e-16",
            "1.96410760376676313e-34",
            "1.99474170985283706e-34",
            "0.165349392879153928",
            "7.12319811999282302e-35",
            "1.45256949579025977e-34",
            "7.91729765076522925e-35",
            "1.40268567204113874e-34",
        ),
        (0, 16, 17, 2, 3, 14),
    )

    # 在同一数值积分器中各取半个电磁读数，第三项保留实际色读数。
    a1, a2, a3 = data["initial"]
    half = (5 * a1 / 3 + a2) / 2
    projected = dict(data, initial=(3 * half / 5, half, a3))
    fields = [(k, z, (q, arb(q) / 2, c)) for k, z, (q, m, c) in species]
    path = lib("fourth-third-pivot-neutral-positive-path")["path"]
    _, floors = path(projected, fields, weights, sites, box, arb(0), lib)
    print("全程电磁、色下界", (2 * floors[0]).str(12), floors[2].str(12), flush=True)
    running = lib("full-neutral-running")
    gauge = lib("all-field-signed-crossing")
    q, m, _ = running["direct"](
        data, sites, weights, species, "3200", gauge, gauge["angular_rule"](13, 39)
    )
    norm = sum(abs(arb(q) - m) for _, _, (q, m, _) in species) / (4 * arb.pi())
    error = norm * sum(
        abs(c) * abs((b / x).log()) for c, b, x in zip(weights, box, sites, strict=True)
    )
    upper = (q - m + weak + error).upper()
    print("3200 GeV 的中性方向 (1,-1) 上界", upper.str(13), flush=True)
    assert upper < arb("-66.76")
    assert 2 * floors[0] > arb("3.6625") and floors[2] > arb("2.0817")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
