"""213 个原壳两面字按字长几何衰减时，两个末阶份额的全区间上下界。"""
from math import comb
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_poly, ctx


def positive_on_unit_interval(poly):
    degree = poly.degree()
    coeff = [poly[k]/comb(degree, k) for k in range(degree+1)]
    bernstein = [sum(comb(j, k)*coeff[k] for k in range(j+1)) for j in range(degree+1)]
    pending = [(bernstein, 0)]
    leaves = 0
    while pending:
        row, depth = pending.pop()
        if all(v > 0 for v in row):
            leaves += 1
            continue
        assert depth < 50
        left, right = [row[0]], [row[-1]]
        while len(row) > 1:
            row = [(u+v)/2 for u, v in zip(row, row[1:])]
            left.append(row[0])
            right.append(row[-1])
        pending.extend(((left, depth+1), (right[::-1], depth+1)))
    return leaves


def bound(label, amplitudes, witness, lower, upper, left, right):
    denominator = sum(a*a for a in amplitudes)
    numerator = amplitudes[-1]**2
    witness, lower, upper, left, right = map(arb, (witness, lower, upper, left, right))
    assert 0 < left < witness < right < 1
    at_witness = numerator(witness)/denominator(witness)
    assert lower < at_witness < upper
    # 上界多项式严格为正也同时认证整个区间归一分母为正。
    leaves = positive_on_unit_interval(upper*denominator-numerator)
    outside = lower*denominator-numerator
    left_leaves = positive_on_unit_interval(outside(arb_poly([0, left])))
    right_leaves = positive_on_unit_interval(outside(arb_poly([right, 1-right])))
    limits = [numerator(t)/denominator(t) for t in (0, 1)]
    print(label, "正性末段数", leaves, left_leaves, right_leaves,
          "见证点份额", at_witness.str(17), "零端与等权端", *(v.str(17) for v in limits))
    return limits


def verify(root):
    geometry = run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))
    ctx.prec = 384
    table = geometry["enumerate_faces"]()
    native = {z for k, z in table.values() if k == (1, 0)}
    first, second = (7514, -3), (29189, -41)
    assert first in native and second in native
    label = geometry["product"](first, geometry["conjugate"](second))
    assert label == (108469, 109)
    u = geometry["unit"](label)
    moments = [arb_poly([1, *[2*(u**(q*j)).real for j in range(1, 107)]])
               for q in range(1, 107)]
    direct = moments[:-1]
    # 差分常数项精确为零。先除去共同正因子 rho，再求连续延拓。
    difference = []
    for a, b in zip(moments, moments[1:]):
        change = b-a
        assert change[0].is_zero()
        difference.append(arb_poly([change[j] for j in range(1, 107)]))
    zero, equal = bound("直接矩", direct, ".99984767", ".01417644590", ".01417644593", ".99984", ".99986")
    assert (zero-arb(1)/105).contains(0) and arb(".01401747") < equal < arb(".01401748")
    zero, equal = bound("相邻差分", difference, ".99964502", ".01428699162", ".01428699165", ".99964", ".99965")
    assert arb(".00000819535") < zero < arb(".00000819536")
    assert arb(".01406353") < equal < arb(".01406354")


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
