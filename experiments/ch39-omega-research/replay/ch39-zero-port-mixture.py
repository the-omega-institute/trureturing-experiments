"""对全正零点谱的非相干混合，用十九个端口的对角读数界定 11431 的权重。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

single = run_path(str(Path(__file__).with_name("ch39-zero-single-port.py")))


def evaluate(poly, x):
    value = arb(0)
    for coefficient in reversed(poly.coeffs()):
        value = value*x + arb(coefficient)
    return value


def verify():
    ctx.prec = 360
    denom, *numerators = single["polynomials"]()
    representatives = [i for i, poly in enumerate(numerators) if poly not in numerators[:i]]
    nums = [numerators[i] for i in representatives]
    slopes = [p.derivative()*denom-p*denom.derivative() for p in nums]
    def row(n):
        height = acb.zeta_zero(n).imag
        x = -1/(2*(height**2+arb(1)/4))
        return [evaluate(p, x)/evaluate(denom, x) for p in nums]
    left, target, right = [row(n) for n in (11430, 11431, 11432)]
    offset = [r-y for r,y in zip(right,target)]
    chord = [l-r for l,r in zip(left,right)]
    domain = arb(-1).union(arb(0))/200
    assert acb.zeta_zero(1).imag > 10 and evaluate(denom, domain) > 0
    for i, j in combinations(range(len(nums)), 2):
        # p 取 11430 的权重；两活跃坐标的残差等幅，法向与邻点弦垂直。
        if chord[i] > 0 and chord[j] > 0 or chord[i] < 0 and chord[j] < 0:
            sign = -1
        else:
            sign = 1
        weight = (sign*offset[j]-offset[i])/(chord[i]-sign*chord[j])
        if not 0 < weight < 1:
            continue
        errors = [b + weight*m for b,m in zip(offset,chord)]
        distance = abs(errors[i])
        if not all(abs(e) < distance for k,e in enumerate(errors) if k not in (i,j)):
            continue
        assert distance > 0 and abs(errors[j]).overlaps(distance)
        # 比较 dy_j/dy_i 的导数。严格同号使邻点弦在全部其他编号的一侧。
        curvature = slopes[i]*slopes[j].derivative()-slopes[j]*slopes[i].derivative()
        assert 0 not in evaluate(curvature, domain) and 0 not in evaluate(slopes[i], domain)
        dual_gap = abs(chord[j]*offset[i]-chord[i]*offset[j])/(abs(chord[i])+abs(chord[j]))
        assert dual_gap.overlaps(distance)
        print("Active ports", representatives[i], representatives[j],
              "11430 weight", weight.str(25), "optimal uniform gap", distance.str(25))
        return
    raise AssertionError("No exposed nearest point certified")


if __name__ == "__main__":
    verify()
