"""精确复现中央十二次多项式的联合子图概率与其边缘乘积之间的总变差和互信息下界。"""
from fractions import Fraction
from itertools import combinations, product
from math import comb, prod
from pathlib import Path
from runpy import run_path
import json

from flint import acb, acb_poly, arb, ctx


REPLAY = Path(__file__).resolve().parent


def geometry_bounds(models, target):
    """在全部联合顶点上包围十二次系数的极差与常数项的最大值。"""
    lows = [None]*12
    highs = [None]*12
    for up, down in product(models[-4], models[2]):
        polynomial = up*down
        for k in range(12):
            lower, upper = polynomial[k].lower(), polynomial[k].upper()
            if lows[k] is None or lower < lows[k]:
                lows[k] = lower
            if highs[k] is None or upper > highs[k]:
                highs[k] = upper
    assert highs[0] < target[0].real**2 + 1
    return max(hi-lo for lo, hi in zip(lows, highs, strict=True))


def root_bounds(factors, models, target):
    """复验每个理想因子分配的分离余量，同时包围六根间距、距离乘积与模。"""
    cover = json.loads((REPLAY / "ch39-independent-geometric-cover.json").read_text())
    coefficients, dot = factors["coefficients"], factors["dot"]
    vertices = {q: [coefficients(f) for f in group] for q, group in models.items()}
    planes = []
    for charge, normal in cover["planes"]:
        assert charge in models and len(normal) == 6 and all(type(x) is int for x in normal)
        norm = sum(map(abs, normal))
        assert norm > 0
        support = max(dot(normal, v).upper() for v in vertices[charge])
        planes.append((charge, normal, norm, support))

    derivatives = [target]
    for _ in range(6):
        derivatives.append(derivatives[-1].derivative())
    as_arb = factors["enclose_rational"]
    triples = list(combinations(range(6), 3))
    separations, products, moduli = [], [], []
    last = Fraction(0)
    for left, right, witnesses in cover["intervals"]:
        left, right = Fraction(left), Fraction(right)
        assert left == last and left < right <= 1
        last = right
        midpoint, half_width = as_arb((right+left)/2), as_arb((right-left)/2)
        roots, boxes, velocities = factors["root_discs"](derivatives, midpoint, half_width)
        moduli.extend(abs(z).upper() for z in boxes)
        for i, box in enumerate(boxes):
            distances = [abs(box-other).lower() for j, other in enumerate(boxes) if i != j]
            separations.append(min(distances))
            products.append(prod(distance/2 for distance in distances))

        for up, witness in zip(triples, witnesses, strict=True):
            assert type(witness) is int and 0 <= witness < len(planes)
            q, normal, norm, support = planes[witness]
            indices = up if q == -4 else tuple(i for i in range(6) if i not in up)
            value, slope = factors["factor_and_derivative"](indices, roots, boxes, velocities)
            margin = dot(normal, value) - support - half_width*abs(dot(normal, slope))
            assert margin > norm*arb("8.74e-7")
    assert last == 1
    return min(separations), min(products), max(moduli)


def verify():
    shell = run_path(str(REPLAY / "ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(REPLAY / "ch39-native-light-projector-pump.py"))
    selector = run_path(str(REPLAY / "ch39-rational-unique-light-action.py"))
    mixture = run_path(str(REPLAY / "ch39-geometric-selector-mixture.py"))
    factors = run_path(str(REPLAY / "ch39-independent-geometric-factors.py"))
    ctx.prec = 320
    target_roots = [arb(t)/9-1 for t in selector["ROOTS"]]
    assert all(abs(a-b) > 0 for a, b in combinations(target_roots, 2))
    target = acb_poly.from_roots([acb(z) for z in target_roots])
    triple = mixture["geometry"](shell, pump)
    models = {q: [mixture["squared_factors"](triple(q, mask, mode))
                  for mask in range(1, 128) if mask != 64 for mode in ("D", "I")]
              for q in (-4, 2)}
    assert all(len(group) == 252 for group in models.values())
    oscillation = geometry_bounds(models, target)
    separation, distance_product, modulus = root_bounds(factors, models, target)
    assert separation > arb(".01631") and distance_product > arb("1.8068e-6")
    assert modulus < arb("1.465") and oscillation < 1566

    # 配对根移动至多 4ρ，m 阶初等对称函数的变化至多为 4ρ·L。
    padded_modulus = modulus + arb("1e-6")
    lipschitz = max(comb(6, m)*m*padded_modulus**(m-1) for m in range(1, 7))
    gap = arb("8.74e-7")
    rho = (gap / (5*lipschitz)).lower()
    assert 4*rho < separation and 4*rho < arb("1e-6") and 4*rho*lipschitz < gap

    # 两个六次因子在圆并边界上各有模下界 ρK。首一多项式相减只余 0 至 11 次。
    coefficient_radius = (rho*distance_product)**2 / sum((modulus+rho)**k for k in range(12))
    total_variation = coefficient_radius / oscillation
    information = 2*total_variation**2
    assert coefficient_radius > arb("1.35e-32")
    assert total_variation > arb("8.67e-36") and information > arb("1.50e-70")
    print("同一六次因子的根间距、五项距离乘积下界、根模上界",
          separation.str(12), distance_product.str(12), modulus.str(12))
    print("纯配方系数最大极差的上界与根圆盘半径", oscillation.str(12), rho.str(12))
    print("十二次系数距离、总变差、互信息（自然对数）下界",
          coefficient_radius.str(12), total_variation.str(12), information.str(12))


if __name__ == "__main__":
    verify()
