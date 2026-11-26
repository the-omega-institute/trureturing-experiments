"""逐枚允许独立共同弱基，对全部四站法向排除左邻与中央同时实现夸克中心。"""
from pathlib import Path
from runpy import run_path
from flint import arb, ctx, fmpq
import numpy as np
ROOT = Path(__file__).resolve().parent
strip = run_path(str(ROOT / 'ch39-left-invariant-character-window.py'))
pencil = run_path(str(ROOT / 'ch39-triple-geometry-pencil.py'))
centred = run_path(str(ROOT / 'ch39-centered-observable-pencil.py'))['centered']
azimuth = run_path(str(ROOT / 'ch39-integer-azimuth-support.py'))
ray_excluded = azimuth['bound_orientation']
raw = run_path(str(ROOT / 'ch39-free-four-site-native-positive-spectrum.py'))
lib, geom, V, mass, A, B, _ = pencil['setup'](ROOT)

def low(x):
    return fmpq(*np.nextafter(float(x.lower()), -np.inf).as_integer_ratio())

def high(x):
    return fmpq(*np.nextafter(float(x.upper()), np.inf).as_integer_ratio())

def cut(poly, a, b, c):
    """以精确有理数裁剪闭半平面 ax+by≤c，保留边界。"""
    if not poly:
        return []
    val = [a * x + b * y - c for x, y in poly]
    if max(val) <= 0:
        return poly
    if min(val) > 0:
        return []
    P = []
    for i, (v, f) in enumerate(zip(poly, val)):
        j = (i + 1) % len(poly)
        if f <= 0:
            P.append(v)
        if f < 0 < val[j] or val[j] < 0 < f:
            r = f / (f - val[j])
            w = poly[j]
            P.append((v[0] + r * (w[0] - v[0]), v[1] + r * (w[1] - v[1])))
    return P

def coeff(a, b, c):
    return (low(a), low(b), high(c))

def limits(poly, a, b, c):
    vs = [a * arb(x) + b * arb(y) + c for x, y in poly]
    return min((i.lower() for i in vs)).union(max((i.upper() for i in vs)))

def polygon(blocks, chart, t):
    """用 G=uv/(1-uv) 与铅笔另一部门的 T22 作非负坐标，求尺度域的外包多边形。"""
    local = chart
    remote = 1 - chart
    d, s, z, k, w, det = blocks[remote]
    if not w > 0:
        return None
    lo_u = max((v.lower() for d_, s_, z_, k_, w_, det_ in blocks for v in (w_ / s_, -w_ / (1 - s_))))
    lo_v = max((v.lower() for d_, s_, z_, k_, w_, det_ in blocks for v in (-w_ / s_, w_ / (1 - s_))))
    if not (lo_u > 0 and lo_v > 0):
        return None
    if lo_u * lo_v > 1:
        return ([], None)
    ratio = -t * mass[local][2] / mass[remote][2]
    dl, sl, zl, kl, wl, delta_l = blocks[local]
    h0 = zl - ratio * z
    h1 = zl - sl - ratio * (z - s)
    h2 = ratio
    D = abs(k) ** 2
    Dl = abs(kl) ** 2
    caps = [(d_ * (1 - d_) / K).upper() for d_, K in ((d, D), (dl, Dl)) if K > 0]
    if not caps:
        return None
    cap = min(caps).fmpq()
    if cap < 1:
        return ([], None)
    poly = [(fmpq(0), fmpq(0)), (cap - 1, fmpq(0)), (cap - 1, fmpq(1)), (fmpq(0), fmpq(1))]
    one = arb(1)
    laws = [(w * (lo_u - 1 / lo_u) - z + s, one, z - w * lo_u), (-w * (1 / lo_v - lo_v) + z - s, -one, -z - w * lo_v), (D / d, -one, -D / d), (D / (1 - d), one, one - D / (1 - d)), (Dl / dl - h1, -h2, h0 - Dl / dl), (Dl / (1 - dl) + h1, h2, one - h0 - Dl / (1 - dl))]
    for aa, bb, cc in laws:
        poly = cut(poly, *coeff(aa, bb, cc))
        if not poly:
            return ([], None)
    return (poly, (h0, h1, h2))

def box_invariants(poly, blocks, h, chart):
    """包围两个迹、上型半径、下型相对分量、上型极角余弦及下型半径。"""
    remote = 1 - chart
    lo = min((x for x, y in poly))
    hi = max((x for x, y in poly))
    L = 1 + arb(lo).union(arb(hi))
    a = arb(min((y for x, y in poly))).union(arb(max((y for x, y in poly))))
    bounds = limits(poly, h[1], h[2], h[0])
    unit_box = arb(0).union(1)
    if not bounds.overlaps(unit_box):
        return False
    a2 = bounds.intersection(unit_box)
    T = [None, None]
    T[remote] = a
    T[chart] = a2
    zz = []
    tr = []
    rs = []
    trans = []
    for i, (d, s, z, k, w, det) in enumerate(blocks):
        tr.append(d + T[i])
        if not tr[-1] > 0:
            return None
        zl = (d.lower() - T[i].upper()) / (d.lower() + T[i].upper())
        zh = (d.upper() - T[i].lower()) / (d.upper() + T[i].lower())
        zz.append(zl.union(zh))
        trans.append(2 * L.sqrt() * k / tr[-1])
        determinant = d * T[i] - L * abs(k) ** 2
        lower = max(arb(0), determinant.lower())
        upper = min(determinant.upper(), (L * det / s).upper())
        if upper < 0:
            return False
        if lower > upper:
            return False
        r2_from_det = 1 - 4 * lower.union(upper) / tr[-1] ** 2
        r2 = zz[-1] ** 2 + abs(trans[-1]) ** 2
        if not r2.overlaps(r2_from_det):
            return False
        r2 = r2.intersection(r2_from_det)
        if not r2 > 0:
            return None
        rs.append(r2.sqrt())
    x, y = trans
    p = (zz[0] * zz[1] + (x * y.conjugate()).real) / rs[0]
    l2 = abs(zz[0] * y - zz[1] * x) ** 2 + (x * y.conjugate()).imag ** 2
    l = arb(0).union(l2.upper()).sqrt() / rs[0] if not l2 > 0 else l2.sqrt() / rs[0]
    if not l.is_finite():
        return None
    polar = zz[0] / rs[0]
    if zz[0] < 0:
        polar = polar.intersection(-1 / (1 + abs(x) ** 2 / zz[0] ** 2).sqrt())
    return (tr[0], tr[1], rs[0], p, l, polar, rs[1])

def up_rays(node):
    up, down = raw['native'](ROOT, 2306888100085 + node)
    rays = strip['native_rays'](up)
    forms = strip['forms'](down)
    ctx.prec = 240
    prep = []
    for j, cap, ray in rays:
        dist = min(((1 - strip['dot'](ray, r)).lower() for k, _, r in rays if k != j))
        first = strip['cross'](ray, [0, 0, 1])
        norm = strip['dot'](first, first).sqrt()
        first = [v / norm for v in first]
        sec = strip['cross'](ray, first)
        tmp = []
        for n0, n, sup in forms:
            v = strip['dot'](n, ray)
            a = strip['dot'](n, first)
            b = strip['dot'](n, sec)
            n2 = strip['dot'](n, n)
            tmp.append((n0, sup, v, a, b, n2))
        prep.append((j, cap, dist, tmp))
    assert len(prep) == 7
    return prep

def certify_poly(poly, h, blocks, chart, sets, budget):
    stack = [(poly, 0)]
    while stack:
        p, level = stack.pop()
        inv = box_invariants(p, blocks, h, chart)
        if inv is False:
            continue
        if inv:
            cut_radius = inv[2].lower()
            high_polar = inv[5].upper()
            axis_gap = (2 * max(arb(0), 1 + high_polar)).sqrt()

            def safe(row, allowed):
                if row[0] not in allowed or not cut_radius > 0:
                    return False
                de = max(arb(0), 1 - cut_radius ** 2)
                K = 2 / row[2]
                pad = (2 * max(arb(0), 1 - cut_radius)).sqrt()
                if K * de < 1:
                    h = (1 - K * de).sqrt()
                    pad = min(pad, (2 * (K - 1) * de ** 2 / (cut_radius * (cut_radius + h) * (1 + h) * (1 + cut_radius))).sqrt().upper())
                return axis_gap + pad < arb('.374')
            need_left = [row for row in sets[-1] if not safe(row, (1, 16, 18))]
            need_central = [row for row in sets[0] if not safe(row, (7, 10, 16, 18))]
            if all((ray_excluded(row, inv) for row in need_left)) and all((ray_excluded(row, inv) for row in need_central)):
                continue
        if level >= budget:
            return False
        a, b, c = (h[1].mid().fmpq(), h[2].mid().fmpq(), h[0].mid().fmpq())
        scalars = [(fmpq(1), 0, 0), (0, fmpq(1), 0), (a, b, c)]
        vals = [[aa * x + bb * y + cc for x, y in p] for aa, bb, cc in scalars]
        widths = [float(max(q) - min(q)) / max(abs(float((max(q) + min(q)) / 2)), [0.1, 0.0001, 0.0001][i]) for i, q in enumerate(vals)]
        j = int(np.argmax(widths))
        aa, bb, cc = scalars[j]
        q = vals[j]
        mid = (max(q) + min(q)) / 2 - cc
        if max(q) == min(q):
            return False
        dyadic = fmpq(*float(mid).as_integer_ratio())
        if min(q) < dyadic + cc < max(q):
            mid = dyadic
        for form in ((aa, bb, mid), (-aa, -bb, -mid)):
            pieces = cut(p, *form)
            if pieces:
                stack.append((pieces, level + 1))
    return True

def prove_chart(chart, idx, scope, depth_limit, sets, geom_depth):
    """在铅笔区间内先覆盖所有尺度。余下的区间继续二分。"""
    base, other = (A, B) if chart == 0 else (B, A)
    stack = [(*scope, 0)]
    vis = done = 0
    mdepth = 0
    while stack:
        l, r, lev = stack.pop()
        t = l.union(r)
        vis += 1
        mdepth = max(lev, mdepth)
        ll = lib['lines'](base, other, l, r, (idx,), geom)
        calc = pencil['physical_blocks'](ll[idx], chart, t, V, mass, geom, lib, (), complex_cross=True) if idx in ll else None
        if calc:
            cc = centred(lib, geom, V, mass, A, B, chart, idx, l, r)
            if cc:
                merged = []
                for raw, new in zip(calc[0], cc):
                    merged.append(tuple((raw[k].intersection(new[k]) if k in (0, 1, 2, 4) else new[k] for k in range(6))))
                calc = (merged,)
        if calc:
            blocks = calc[0]
            P = polygon(blocks, chart, t)
            if P is not None:
                poly, h = P
                if not poly or certify_poly(poly, h, blocks, chart, sets, geom_depth):
                    done += 1
                    continue
        if lev >= depth_limit:
            print('尚未覆盖的坐标图、本征支和参数区间', chart, idx, vis, float(l), float(r), flush=True)
            return False
        mid = (l + r) / 2
        stack.extend([(mid, r, lev + 1), (l, mid, lev + 1)])
    print('图、支、铅笔盒、成功叶盒、最大深度', chart, idx, vis, done, mdepth, flush=True)
    return True

def verify():
    azimuth['compile_cover'](ROOT)
    sets = {node: up_rays(node) for node in (-1, 0)}
    for chart in (0, 1):
        for branch in range(3):
            if not prove_chart(chart, branch, (arb(-1), arb(1)), 22, sets, 14):
                raise ArithmeticError("尚有未覆盖的铅笔区间")
    left = strip['native_rays'](raw['native'](ROOT, 2306888100084)[0])
    centre = strip['native_rays'](raw['native'](ROOT, 2306888100085)[0])
    minimum = min((sum(((a - b) ** 2 for a, b in zip(x, y, strict=True))).sqrt().lower() for i, _, x in left if i in (1, 16, 18) for j, _, y in centre if j in (7, 10, 16, 18)))
    assert minimum > arb('.7504203386270') > 2 * arb('.374')
    print('两组实际射线最小弦长、共同轴允许的两倍半径', minimum.str(16), 2 * arb('.374'))
if __name__ == '__main__':
    verify()
