"""两荷独立抽取无钟子图时，期望多项式不能精确等于校准多项式加实常数。

取 x=t/9−1，p 为六个校准根对应的首一六次多项式。独立性将期望分解为
两个首一实六次多项式的乘积，每个因子在本荷 252 种三根平方积的凸包中。
C>0 时，p²+C 的实六次因子必须从六对共轭根中选三对。C=0 时，两个
期望因子都处处非负，p 的每个单根在乘积中的重数二必须全部归于一方。
C<0 时，在 p 的任一实根处，p²+C 为负。

同目录 JSON 的 planes 各行给出荷与六个整数系数。intervals 各行给出
s=√C 的有理端点与二十个从零起算的面编号。分法依 combinations 顺序，
每次列出的三对归上型，补集归下型。法向量按 x⁵ 至常数项的顺序作用。
"""
from fractions import Fraction
from itertools import combinations
from math import factorial
from pathlib import Path
from runpy import run_path
import json

from flint import acb, acb_poly, arb, ctx


REPLAY = Path(__file__).resolve().parent
TRIPLES = tuple(combinations(range(6), 3))
ALL_ROOTS = set(range(6))


def enclose_rational(value):
    return arb(value.numerator) / value.denominator


def coefficients(polynomial):
    return [polynomial[k].real for k in range(5, -1, -1)]


def conjugate(polynomial):
    return acb_poly([z.conjugate() for z in polynomial.coeffs()])


def dot(left, right):
    return sum(x*y for x, y in zip(left, right, strict=True))


def subgraph_coefficients(triple, squared_factors, charge):
    coeffs = []
    # 六个三站轨道占低六位，极点占第七位。
    for mask in range(1, 128):
        if mask == 64:
            continue
        for mode in ("D", "I"):
            roots = triple(charge, mask, mode)
            coeffs.append(coefficients(squared_factors(roots)))
    assert len(coeffs) == 252
    return coeffs


def root_discs(derivatives, midpoint, half_width):
    """由 Rouché 定理在整个 s 区间连续跟踪六个互不相交的根。"""
    middle_polynomial = derivatives[0] + acb(0, midpoint)
    roots = sorted(middle_polynomial.roots(tol=arb(2)**-180, maxprec=1024),
                   key=lambda z: float(z.real.mid()))
    assert len(roots) == 6
    centers, radii, boxes = [], [], []
    for root in roots:
        center = acb(root.real.mid(), root.imag.mid())
        residual = abs(middle_polynomial(center)).upper() + half_width
        slope = abs(derivatives[1](center)).lower()
        assert slope > 0
        radius = (arb("1.4") * residual / slope).upper()
        tail = sum(abs(derivatives[k](center)) * radius**k / factorial(k)
                   for k in range(2, 7))
        assert abs(root-center) < radius and slope*radius > residual+tail
        centers.append(center)
        radii.append(radius)
        boxes.append(center + acb(arb(0, radius), arb(0, radius)))

    for i, j in combinations(range(6), 2):
        assert abs(centers[i]-centers[j]) > radii[i]+radii[j]
    # p(z)+is=0 给出 z′=−i/p′(z)。以圆盘外接矩形包围导数。
    velocities = [-acb(0, 1)/derivatives[1](box) for box in boxes]
    assert all(z.is_finite() for z in velocities)
    return roots, boxes, velocities


def factor_and_derivative(indices, roots, boxes, velocities):
    center = acb_poly.from_roots([roots[j] for j in indices])
    middle_coeffs = coefficients(center * conjugate(center))
    factors = [acb_poly([-boxes[j], 1]) for j in indices]
    product = factors[0]*factors[1]*factors[2]
    derivative = -(velocities[indices[0]] * factors[1]*factors[2]
                   + velocities[indices[1]] * factors[0]*factors[2]
                   + velocities[indices[2]] * factors[0]*factors[1])
    derivative_coeffs = coefficients(derivative*conjugate(product) + product*conjugate(derivative))
    return middle_coeffs, derivative_coeffs


def verify():
    shell = run_path(str(REPLAY / "ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(REPLAY / "ch39-native-light-projector-pump.py"))
    mixture = run_path(str(REPLAY / "ch39-geometric-selector-mixture.py"))
    selector = run_path(str(REPLAY / "ch39-rational-unique-light-action.py"))
    ctx.prec = 320

    triple = mixture["geometry"](shell, pump)
    models = {q: subgraph_coefficients(triple, mixture["squared_factors"], q) for q in (-4, 2)}
    target_roots = [arb(r)/9-1 for r in selector["ROOTS"]]
    assert all(abs(a-b) > 0 for a, b in combinations(target_roots, 2))
    p = acb_poly.from_roots([acb(r) for r in target_roots])
    # C≥1 的常数项超过两个凸组合常数项的乘积上界。各顶点均为实多项式平方。
    maxima = [max(v[-1].upper() for v in models[q]) for q in (-4, 2)]
    assert p[0].real**2 + 1 > maxima[0]*maxima[1]
    derivatives = [p]
    for _ in range(6):
        derivatives.append(derivatives[-1].derivative())

    data = json.loads((REPLAY / "ch39-independent-geometric-cover.json").read_text())
    planes = []
    for charge, normal in data["planes"]:
        assert charge in models and len(normal) == 6 and all(type(x) is int for x in normal)
        norm = sum(abs(x) for x in normal)
        assert norm > 0
        support = max(dot(normal, v).upper() for v in models[charge])
        planes.append((charge, normal, norm, support))

    last = Fraction(0)
    minimum = None
    count = 0
    for left, right, witnesses in data["intervals"]:
        left, right = Fraction(left), Fraction(right)
        assert left == last and left < right <= 1
        last = right
        midpoint = enclose_rational((left+right)/2)
        half_width = enclose_rational((right-left)/2)
        roots, boxes, velocities = root_discs(derivatives, midpoint, half_width)
        for up, witness in zip(TRIPLES, witnesses, strict=True):
            assert type(witness) is int and 0 <= witness < len(planes)
            charge, normal, norm, support = planes[witness]
            indices = up if charge == -4 else tuple(sorted(ALL_ROOTS.difference(up)))
            value, derivative = factor_and_derivative(indices, roots, boxes, velocities)
            # 中点的分离余量减去整段上函数变化的上界，余量按法向量一范数归一。
            gap = (dot(normal, value) - support - half_width*abs(dot(normal, derivative))).lower()/norm
            assert gap > 0
            minimum = gap if minimum is None else min(minimum, gap)
            count += 1
    assert last == 1 and len(TRIPLES) == 20 and count == 20*len(data["intervals"])
    print("连续闭区间、已验证分法、整数分离面", len(data["intervals"]), count, len(planes))
    print("一范数归一后的最小严格分离余量", minimum.str(12))


if __name__ == "__main__":
    verify()
