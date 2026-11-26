"""同一端口相位参考下，整个十九点逆响应区分 11431 号零点的锐界。"""
from collections import Counter
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

shell = run_path(str(Path(__file__).with_name("ch39-refined-mass-obstruction.py")))


def verify():
    ctx.prec = 192
    _, edges = shell["geometry"]()
    charged = [(i, j, w, q) for i, j, w, q in edges if q]
    weight = Fraction(27, 52136)
    assert all(j < 19 and w == weight and abs(q) == 1 for i, j, w, q in charged)
    neighbors = {}
    for i, j, _, _ in charged:
        neighbors.setdefault(i, set()).add(j)
        neighbors.setdefault(j, set()).add(i)
    centers = [i for i, row in neighbors.items() if len(row) == 2]
    assert Counter(map(len, neighbors.values())) == {1: 6, 2: 3}
    assert all(len(neighbors[j]) == 1 for i in centers for j in neighbors[i])
    heights = []
    phases = []
    for n in (11430, 11431, 11432):
        rho = acb.zeta_zero(n)
        assert rho.real == arb(1)/2
        heights.append(rho.imag)
        phases.append((rho-1)/rho)
    radii = []
    for i in (0, 2):
        chord = abs(phases[1]-phases[i])
        a, b = heights[1], heights[i]
        analytic = abs(a-b)/((a*a+arb(1)/4)*(b*b+arb(1)/4)).sqrt()
        assert analytic.overlaps(chord)
        radii.append(arb(2).sqrt()*shell["ball"](weight)*chord/2)
    left, right = radii
    assert arb("3.33210361251700e-12") < left < arb("3.33210361251702e-12")
    assert arb("1.79399759096856e-12") < right < arb("1.79399759096859e-12")
    print("Target 11431, direct inverse-response operator norm, left/right half gaps:",
          left.str(22), right.str(22), flush=True)


if __name__ == "__main__":
    verify()
