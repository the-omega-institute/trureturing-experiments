"""五个连续编号中的三个两两不相容的夸克空间要求。"""
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import arb
S = 10000000
LAYERS = {-2: (1, 7, 10, 16, 18), -1: (16, 18), 0: (7, 10, 16, 18)}
CHORD = arb('.225')

def lower(x):
    return int((S * x).lower().floor().unique_fmpz())

def upper(x):
    return int((S * x).upper().ceil().unique_fmpz())

def bounds(x):
    return (lower(x), upper(x))

def product(a, b, c, d):
    t = (a * c, a * d, b * c, b * d)
    return (np.minimum.reduce(t), np.maximum.reduce(t))

def axis_components(model, poly, h, blocks, chart, inv):
    """在上型轴、下型横向和有向叉积的正交基中包围物理轴 −Z 的三分量。"""
    tu, td, r, P, l, polar, rd = inv
    if not l > 0:
        return None
    remote = 1 - chart
    q = 1 + arb(min((x for x, y in poly))).union(arb(max((x for x, y in poly))))
    rem = arb(min((y for x, y in poly))).union(arb(max((y for x, y in poly))))
    loc = model['limits'](poly, h[1], h[2], h[0]).intersection(arb(0).union(1))
    z = [None, None]
    z[chart] = loc
    z[remote] = rem
    u, d = blocks
    du, dd = (u[0], d[0])
    ku, kd = (u[3], d[3])
    twice = 4 * q * (ku * kd.conjugate())
    h1 = ((du - z[0]) * twice.real - (dd - z[1]) * 4 * q * abs(ku) ** 2) / (tu ** 2 * td * r ** 2 * l)
    h2 = -twice.imag / (tu * td * r * l)
    physical = arb(-1).union(arb(1))
    result = []
    for w in (-polar, h1, h2):
        if not w.is_finite():
            return None
        if not w.overlaps(physical):
            return False
        result.append(w.intersection(physical))
    return result

def prepare(model, node):
    strip = model['strip']
    up, down = model['raw']['native'](model['ROOT'], 2306888100085 + node)
    rays = strip['native_rays'](up)
    supports = strip['forms'](down)
    frames = {}
    rows = []
    for j, capacity, ray in rays:
        near = min(((1 - strip['dot'](ray, r)).lower() for k, _, r in rays if k != j))
        f = strip['cross'](ray, [0, 0, 1])
        length = strip['dot'](f, f).sqrt()
        f = [x / length for x in f]
        g = strip['cross'](ray, f)
        forms = [(n0, top, strip['dot'](n, ray), strip['dot'](n, f), strip['dot'](n, g), strip['dot'](n, n)) for n0, n, top in supports]
        frames[j] = (ray[0], f[0], g[0], {})
        rows.append((j, capacity, near, forms))
    assert len(rows) == 7
    return (rays, rows, frames)

def in_cap(row, radius, polar):
    if not radius > 0:
        return False
    defect = max(arb(0), 1 - radius ** 2)
    K = 2 / row[2]
    pad = (2 * max(arb(0), 1 - radius)).sqrt()
    if K * defect < 1:
        h = (1 - K * defect).sqrt()
        pad = min(pad, (2 * (K - 1) * defect ** 2 / (radius * (radius + h) * (1 + h) * (1 + radius))).sqrt().upper())
    return (2 * max(arb(0), 1 + polar)).sqrt() + pad < CHORD

def positive_segment(model, frame, depth):
    ex, fx, gx, cache = frame
    az = model['azimuth']
    if depth not in cache:
        a = []
        b = []
        for j in range(2 ** depth):
            sn, cs = az['circle'](depth, j)
            a.append(bounds(fx * cs + gx * sn))
            b.append(bounds(gx * cs - fx * sn))
        cache[depth] = (np.array(a, np.int64).T, np.array(b, np.int64).T)
    return cache[depth]

def axis_or_violation(model, row, frame, inv, components, side):
    """对每段角度证明所需的有向坐标界，或证明下型原生支撑严格违约。"""
    az = model['azimuth']
    j, cap, dist, forms = row
    tu, td, r, p, l, polar, rd = inv
    de = 1 - r ** 2
    if dist * (1 - cap / tu) > de:
        return True
    fac = 2 / dist
    if not (r > 0 and fac * de < 1):
        return False
    h = (1 - fac * de).sqrt()
    loss = ((fac - 1) * de ** 2 / (r * (r + h) * (1 + h) * (1 + r))).upper()
    pad = (2 * loss).sqrt()
    if not all((x.is_finite() for x in (r, p, l, rd, pad))):
        return False
    if components is False:
        return True
    if components is None:
        return False
    scaled = tuple((bounds(x) for x in components))
    q = (lower(p), upper(p), lower(l), upper(l), lower(td), upper(rd), upper(loss), upper(pad))
    if q[4] <= 0 or max(map(abs, q)) > 4 * S:
        return False
    E, F, G = scaled
    am, ap = product(*E, *bounds(frame[0]))
    limit = (-3000000 if side > 0 else -3200000) * S
    gap = upper(pad) * S
    flat, wide = az['initialize'](forms)
    pointer = az['point']
    integer = az['compile_cover'].__globals__['support_library']
    box = np.array(q, np.int64)
    for depth in (4, 6, 8):
        a, b = positive_segment(model, frame, depth)
        lm, lp = product(*F, *a)
        mm, mp = product(*G, *b)
        ok = am + lm + mm - gap > limit if side > 0 else ap + lp + mp + gap < limit
        if np.all(ok):
            return True
        spans = np.ascontiguousarray(az['section'](forms, depth, wide)[~ok])
        if integer.supports(len(forms), len(spans), pointer(flat), pointer(spans), pointer(box)):
            return True
    return False

def cover_scales(model, outer, initial, diag, blocks, chart, spectra, frames):
    """在有理尺度多边形上核对每个可能层的弦长和两组所需的有向坐标界。"""
    stack = [(initial, 0)]
    while stack:
        poly, depth = stack.pop()
        inv = model['box_invariants'](poly, blocks, diag, chart)
        if inv is False:
            continue
        if inv:
            radius, polar = (inv[2].lower(), inv[5].upper())
            passed = True
            components = axis_components(model, poly, diag, blocks, chart, inv)
            for node, rows in spectra.items():
                for row in rows:
                    j = row[0]
                    side = 1 if node == -2 and j in (7, 10) else -1 if node == 0 and j == 7 else 0
                    if j in LAYERS[node] and in_cap(row, radius, polar):
                        if not side or axis_or_violation(model, row, frames[node][j], inv, components, side):
                            continue
                    if model['ray_excluded'](row, inv):
                        continue
                    passed = False
                    break
                if not passed:
                    break
            if passed:
                continue
        if depth >= 14:
            return False
        pieces = outer['split_polygon'](model, poly, diag)
        if not pieces:
            return False
        stack.extend(((v, depth + 1) for v in pieces if v))
    return True

def cover_branch(model, outer, chart, branch, spectra, frames):
    pencil, lib, geom = (model['pencil'], model['lib'], model['geom'])
    mixing, mass, a, b = (model[key] for key in ('V', 'mass', 'A', 'B'))
    base, other = (a, b) if chart == 0 else (b, a)
    stack = [(arb(-1), arb(1), 0)]
    visits = leaves = deepest = 0
    while stack:
        low, high, depth = stack.pop()
        visits += 1
        deepest = max(depth, deepest)
        interval = low.union(high)
        lines = lib['lines'](base, other, low, high, (branch,), geom)
        if branch in lines:
            calc = pencil['physical_blocks'](lines[branch], chart, interval, mixing, mass, geom, lib, (), complex_cross=True)
            if calc:
                blocks = calc[0]
                refined = model['centred'](lib, geom, mixing, mass, a, b, chart, branch, low, high)
                if refined:
                    blocks = [tuple((old[k].intersection(new[k]) if k in (0, 1, 2, 4) else new[k] for k in range(6))) for old, new in zip(blocks, refined, strict=True)]
                result = model['polygon'](blocks, chart, interval)
                if result is not None:
                    poly, diagonal = result
                    if not poly or cover_scales(model, outer, poly, diagonal, blocks, chart, spectra, frames):
                        leaves += 1
                        continue
        if depth >= 22:
            raise ArithmeticError((chart, branch, low, high))
        middle = (low + high) / 2
        stack.extend(((middle, high, depth + 1), (low, middle, depth + 1)))
    print('图、支、铅笔盒、成功叶盒、最大深度', chart, branch, visits, leaves, deepest, flush=True)

def verify(root):
    model = run_path(str(root / 'ch39-global-independent-weak-frame-separation.py'))
    outer = run_path(str(root / 'ch39-left-right-independent-weak-frame-separation.py'))
    model['azimuth']['compile_cover'](root)
    rays, rows, frames = ({}, {}, {})
    for node in LAYERS:
        rays[node], rows[node], frames[node] = prepare(model, node)
    for chart in (0, 1):
        for branch in range(3):
            cover_branch(model, outer, chart, branch, rows, frames)
    exceptional = {(7, 7), (10, 7)}
    for i, j in ((-2, -1), (-1, 0), (-2, 0)):
        excluded = exceptional if (i, j) == (-2, 0) else set()
        distances = {(k, l): sum(((a - b) ** 2 for a, b in zip(x, y, strict=True))).sqrt() for k, _, x in rays[i] if k in LAYERS[i] for l, _, y in rays[j] if l in LAYERS[j]}
        smallest = min((d.lower() for key, d in distances.items() if key not in excluded))
        assert smallest > arb('.450')
        if excluded:
            assert {key for key, d in distances.items() if d < arb('.450')} == excluded
        print('两个壳偏移、弦长直接排除的最小距离', i, j, smallest.str(16))
    assert arb('-.30') > arb('-.32')
if __name__ == '__main__':
    verify(Path(__file__).resolve().parent)
