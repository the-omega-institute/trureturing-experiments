"""用两枚外邻各自的真实谱层检验可行物理轴是否相交。"""
from pathlib import Path
from runpy import run_path

from flint import arb, fmpq


ROOT = Path(__file__).resolve().parent
CENTRAL = 2_306_888_100_085
RADIUS = ".335"
LAYERS = {-1: (16, 18), 1: (1, 7, 10, 16, 18)}


def small_cap(row, radius, polar):
    """若本层为最近射线，界定该射线到物理方向 -Z 的弦长。"""
    if not radius > 0:
        return False
    defect = max(arb(0), 1 - radius**2)
    factor = 2 / row[2]
    pad = (2 * max(arb(0), 1-radius)).sqrt()
    if factor * defect < 1:
        root = (1-factor*defect).sqrt()
        loss = (factor-1)*defect**2 / (radius*(radius+root)*(1+root)*(1+radius))
        pad = min(pad, (2*loss).sqrt().upper())
    return (2 * max(arb(0), 1+polar)).sqrt() + pad < arb(RADIUS)


def split_polygon(model, poly, diagonal):
    """按 G 或两块的第二对角元的相对宽度，选一条有理割线。"""
    c, a, b = (value.mid().fmpq() for value in diagonal)
    coordinates = ((1, 0, 0), (0, 1, 0), (a, b, c))
    ranges = []
    for (ax, ay, start), scale in zip(coordinates, (.1, .0001, .0001), strict=True):
        values = [ax*x + ay*y + start for x, y in poly]
        low, high = min(values), max(values)
        width = float(high-low) / max(scale, abs(float((high+low)/2)))
        ranges.append((width, low-start, high-start, ax, ay))
    width, low, high, ax, ay = max(ranges, key=lambda entry: entry[0])
    if not width:
        return ()
    middle = (low+high)/2
    dyadic = fmpq(*float(middle).as_integer_ratio())
    if low < dyadic < high:
        middle = dyadic
    return (model["cut"](poly, ax, ay, middle), model["cut"](poly, -ax, -ay, -middle))


def cover_scales(model, initial, diagonal, blocks, chart, spectra):
    stack = [(initial, 0)]
    while stack:
        poly, depth = stack.pop()
        inv = model["box_invariants"](poly, blocks, diagonal, chart)
        if inv is False:
            continue
        if inv:
            radius, polar = inv[2].lower(), inv[5].upper()
            passes = all(
                (row[0] in LAYERS[offset] and small_cap(row, radius, polar))
                or model["ray_excluded"](row, inv)
                for offset, rows in spectra.items() for row in rows
            )
            if passes:
                continue
        if depth >= 14:
            return False
        pieces = split_polygon(model, poly, diagonal)
        if not pieces:
            return False
        stack.extend((piece, depth+1) for piece in pieces if piece)
    return True


def cover_branch(model, chart, branch, spectra):
    pencil, lib, geom = model["pencil"], model["lib"], model["geom"]
    mixing, mass, a, b = (model[key] for key in ("V", "mass", "A", "B"))
    base, other = (a, b) if chart == 0 else (b, a)
    stack = [(arb(-1), arb(1), 0)]
    visits = leaves = deepest = 0
    while stack:
        low, high, depth = stack.pop()
        visits += 1
        deepest = max(depth, deepest)
        interval = low.union(high)
        lines = lib["lines"](base, other, low, high, (branch,), geom)
        if branch in lines:
            calc = pencil["physical_blocks"](lines[branch], chart, interval, mixing, mass,
                                            geom, lib, (), complex_cross=True)
            if calc:
                blocks = calc[0]
                refined = model["centred"](lib, geom, mixing, mass, a, b, chart, branch, low, high)
                if refined:
                    blocks = [tuple(old[k].intersection(new[k]) if k in (0, 1, 2, 4)
                                    else new[k] for k in range(6))
                              for old, new in zip(blocks, refined, strict=True)]
                result = model["polygon"](blocks, chart, interval)
                if result is not None:
                    poly, diagonal = result
                    if not poly or cover_scales(model, poly, diagonal, blocks, chart, spectra):
                        leaves += 1
                        continue
        if depth >= 22:
            raise ArithmeticError((chart, branch, low, high))
        middle = (low+high)/2
        stack.extend(((middle, high, depth+1), (low, middle, depth+1)))
    print("图、支、铅笔盒、成功叶盒、最大深度", chart, branch, visits, leaves, deepest, flush=True)


def verify(root):
    model = run_path(str(root/"ch39-global-independent-weak-frame-separation.py"))
    model["azimuth"]["compile_cover"](root)
    spectra = {offset: model["up_rays"](offset) for offset in LAYERS}
    for chart in (0, 1):
        for branch in range(3):
            cover_branch(model, chart, branch, spectra)
    strip, native = model["strip"], model["raw"]["native"]
    sets = {}
    for offset, allowed in LAYERS.items():
        rows = strip["native_rays"](native(root, CENTRAL+offset)[0])
        sets[offset] = [(i, ray) for i, _, ray in rows if i in allowed]
    minimum, i, j = min((sum((a-b)**2 for a, b in zip(x, y, strict=True)).sqrt().lower(), i, j)
                        for i, x in sets[-1] for j, y in sets[1])
    assert minimum > arb(".6757960768330") > 2*arb(RADIUS)
    print("最近左、右原生层、最小弦长、同一物理轴的弦长上界", i, j, minimum.str(16), 2*arb(RADIUS))


if __name__ == "__main__":
    verify(ROOT)
