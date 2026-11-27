"""第五壳以三个高截止和十六个低截止同时实现五项条件与全程正动能。"""
from pathlib import Path
from runpy import run_path

from flint import arb


CENTERS = (
    "7.45405576985424e-35", "8.93645866900793e-23", ".08180490720755705",
    ".2424980115609757", "8.927421430262749e-23", "8.716300564362518e-23",
    "8.707517119074346e-23", "8.896034056853818e-23", "9.063339869188606e-23",
    "6.865414589285544e-23", "6.00915444274547e-23", "8.858037012788777e-23",
    "9.006513620484771e-23", "8.936300520925373e-23", ".16534944339482094",
    "5.107092418604431e-34", "1.0432114428059044e-34", "2.5829186554312497e-34",
    "3.451030012462339e-30",
)
AXES = (17, 0, 15, 2, 3, 14)


def endpoints(root):
    solve = run_path(str(root / "ch39-fifth-third-pivot-em-colour-path.py"))["solve"]
    data, species, weights, centers, boxes, _, lib = solve(root, CENTERS, AXES)
    assert [j+1 for j, x in enumerate(boxes) if x > arb("1e-6")] == [3, 4, 15]
    assert all(x < arb("1e-6") for j, x in enumerate(boxes) if j not in (2, 3, 14))
    return data, species, weights, centers, boxes, lib


def weak_along(data, weights, boxes):
    ratio = (arb("246.22") / arb("91.1876"))**2 / 2
    w = arb(".63811")**2 * ratio / 2
    z = w + 3 * arb(".467453")**2 * ratio / 10
    t = arb(".8616")**2 * ratio
    boundaries = [data["p0"] / x for x in boxes]
    lo, hi = arb(2)**-14, arb(8)
    # 在短固有时段用 |C1|+C2 的一致界，长时段用已归一的尾轮廓。
    early = sum(map(abs, weights)) * lo * (5*w + 3*z + 5*t) / (48*arb.pi())
    late = (4 + 9 / arb(1).exp()) / (8 * arb.pi() * hi)
    assert all(q < hi for q in boundaries) and sum(weights).contains(1)
    nodes = sorted([arb(2)**j for j in range(-14, 4)] +
                   [v for v in boundaries if v > lo], key=lambda v: v.mid())
    assert nodes[0] == lo and nodes[-1] == hi
    entries = []
    angles = 32
    for a, b in zip(nodes, nodes[1:]):
        assert a < b
        middle = (a + b) / 2
        assert all(v < middle or v > middle for v in boundaries)
        profile = sum(c for c, v in zip(weights, boundaries, strict=True) if v < middle)
        steps = 6
        for k in range(steps):
            cell = (a + (b-a)*k/steps).union(a + (b-a)*(k+1)/steps)
            common = 4*w*(-w*cell).exp() - 3*t*(-t*cell).exp()
            # (1-exp(-t u))/u 单调递减。
            vector = (-2*(-t*cell.lower()).expm1()/cell.lower()).union(
                      -2*(-t*cell.upper()).expm1()/cell.upper())
            mass_factor = profile * (b-a)/(steps*angles*8*arb.pi())
            for i in range(angles):
                left, right = arb(i)/angles, arb(i+1)/angles
                sh = left.union(right)
                h = (left*(1-left)).union(right*(1-right))
                angular_mass = sh*z
                # (exp(-s z u)-exp(-w u))/u 的平均值定理。
                neutral = (w-angular_mass)*(-cell*w).exp().union((-cell*angular_mass).exp())
                scalar = common + 2*z*(-angular_mass*cell).exp() + neutral
                y = h*cell
                entries.append((max(arb(0), y.lower()), y.upper(),
                                mass_factor*h*scalar, mass_factor*h*vector))

    def values(y):
        if y.is_zero():
            return arb(1), arb(1)
        return -(-y).expm1()/y, (-y).exp()

    # Φ(P y) 与 exp(-P y) 对 y 与 P 都递减，逐矩形取两个对角。
    def evaluate(pl, ph):
        result = arb(0)
        for yl, yh, first, second in entries:
            f_hi, e_hi = values(pl*yl)
            if ph is None:
                f_lo = e_lo = arb(0)
            else:
                f_lo, e_lo = values(ph*yh)
            result += first * f_lo.union(f_hi) + second * e_lo.union(e_hi)
        return result + arb(0, (early+late/pl).upper())

    # 对 P1 的独立包围也用于同一个精确根的倒推。
    end = data["p1"]/data["p0"]
    fixed = evaluate(end, end)
    span = None
    for k in range(12):
        for j in range(4):
            pl = arb(2)**k*(4+j)/4
            ph = arb(2)**k*(5+j)/4
            bound = evaluate(pl, ph) - fixed
            span = bound if span is None else span.union(bound)
    bound = evaluate(arb(2)**12, None) - fixed
    span = span.union(bound)
    print("全程弱差包围、无穷参考剩余", span.lower().str(12), span.upper().str(12),
          fixed.str(12), flush=True)
    assert arb("-.448876") < span < arb(".732559")
    return -span.lower()


def verify(root):
    data, species, weights, centers, boxes, lib = endpoints(root)
    weak = weak_along(data, weights, boxes)
    # 色读数在调用共享程序时乘二，得到正的物理色下界后再除二。
    a1, a2, a3 = data["initial"]
    scaled = dict(data, initial=(a1, a2, 2*a3))
    rows = [(k,z,(q,m,2*c)) for k,z,(q,m,c) in species]
    path = lib("fourth-third-pivot-neutral-positive-path")["path"]
    _, floors = path(scaled, rows, weights, centers, boxes, weak, lib)
    assert floors[2]/2 > arb("1.00938") and floors[3] > arb("8.29830")
    print("连续色逆耦合、行列式下界", (floors[2]/2).str(13), floors[3].str(13))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
