"""用十九点球壳磁拉普拉斯的有限行列式生成候选面权，检验三支迹比。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_poly, ctx


def determinant_polynomial(matrix):
    coefficients = list((-matrix).charpoly())
    assert all(value.imag.contains(0) for value in coefficients)
    return arb_poly([value.real for value in coefficients])


def stiffness(polynomials, x, family, dimension, index):
    totals = [arb(0) for _ in range(3)]
    for rep in family:
        p = polynomials[abs(rep[0])]
        multiplier = 2 * p.derivative()(x) / p(x)
        for j, group_index in enumerate(index(rep)):
            totals[j] += multiplier * dimension(rep) * group_index.numerator / group_index.denominator
    return totals


def signed_numerator(polynomials, weights, sign):
    # P'_q / P_q = 19/x - tr(L_q)/x² + tr(L_q²)/x³ + O(x^-4).
    # The first two traces depend only on edge magnitudes. The weighted
    # combination cancels all three displayed terms exactly.
    assert sum(weights.values()) == 0
    numerator = sum(a * polynomials[q].derivative()
                    * prod(polynomials[k] for k in weights if k != q)
                    for q, a in weights.items())
    degree = 19 * len(weights) - 4
    assert all(numerator[j].contains(0) for j in range(degree + 1, numerator.degree() + 1))
    assert all(sign * numerator[j] > 0 for j in range(degree + 1))
    return degree


def verify():
    directory = Path(__file__).parent
    shell = run_path(str(directory / "ch39-three-ray-chiral-boundary.py"))
    trace = run_path(str(directory / "ch39-native-trace-plaquette-jets.py"))
    family, dimension, index = (trace[k] for k in ("family", "dimension", "index"))
    ctx.prec = 512
    height = acb.zeta_zero(11431).imag
    charges = [abs(rep[0]) for rep in family]
    polynomials = {q: determinant_polynomial(shell["laplacian"](q, height)) for q in charges}
    assert all(p.degree() == 19 and all(a > 0 for a in p) for p in polynomials.values())
    flat = determinant_polynomial(shell["laplacian"](0, height))
    flat_y, flat_w, flat_c = stiffness(dict.fromkeys(charges, flat), 4, family, dimension, index)
    assert (flat_w - flat_c).contains(0) and (flat_y - 60 * flat_w).contains(0)
    y, w, c = stiffness(polynomials, 4, family, dimension, index)
    assert arb(".99561348385") < w / c < arb(".99561348387")
    assert arb("1.01493830180") < y / (60 * w) < arb("1.01493830182")
    print("弱/色、超荷/(60 弱)", (w / c).str(14), (y / (60 * w)).str(14))
    # These two combinations have the signs of W-C and Y-60W, respectively.
    assert signed_numerator(polynomials, {1: 1, 3: 1, 2: -1, 4: -1}, -1) == 72
    assert signed_numerator(polynomials, {1: -84, 2: 12, 3: -12, 4: 48, 6: 36}, 1) == 91
    theta = shell["geo"]["theta"]
    for d, lower, upper in ((724, "1.96215e-24", "1.96216e-24"),
                            (804, "1.46434e-24", "1.46435e-24")):
        angle = height * (theta(d) - theta(244))
        face = arb(1)
        for rep in family:
            q = rep[0]
            p = polynomials[abs(q)]
            face *= (p(2 + 2 * (q * angle).cos()) / p(4))**dimension(rep)
        assert arb(lower) < face < arb(upper)
        print("第 11431 枚零点、弦类与候选超荷面权", d, face.str(14))
    top_frequency = sum(19 * abs(rep[0]) * dimension(rep) for rep in family)
    assert top_frequency == 684


if __name__ == "__main__":
    verify()
