"""逐个比较十九个旧点在时钟 1 的单端口对角预解读数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpq, fmpq_mat, fmpq_poly

memory = run_path(str(Path(__file__).with_name("ch39-three-clock-minimal-memory.py")))


def polynomials():
    _, _, _, g0, g1, g2 = memory["coefficients"]()
    _, edges = memory["REFINE"]["geometry"]()
    base = memory["identity"](19)-g0-g1-g2
    charged = []
    for i, j, weight, exponent in edges:
        if j >= 19:
            continue
        weight = fmpq(weight.numerator, weight.denominator)
        if exponent:
            charged.append((i, j, weight, exponent))
        else:
            base[i, j] -= weight
            base[j, i] -= weight
    assert len(charged) == 6
    xs, values = [], []
    for sample in range(1, 8):
        w = fmpq(sample)
        xs.append((w+1/w)/2-1)
        pencil = fmpq_mat(base)
        for i, j, weight, exponent in charged:
            pencil[i, j] -= weight*w**exponent
            pencil[j, i] -= weight*w**(-exponent)
        determinant = pencil.det()
        inverse = pencil.inv()
        values.append([determinant, *(determinant*inverse[i, i] for i in range(19))])
    polys = [fmpq_poly() for _ in values[0]]
    for i, row in enumerate(values):
        basis = fmpq_poly([1])
        for j, x in enumerate(xs):
            if i != j:
                basis *= fmpq_poly([-x, 1])/(xs[i]-x)
        for k, value in enumerate(row):
            polys[k] += basis*value
    return polys


def verify():
    ctx.prec = 240
    denom, *nums = polynomials()
    assert acb.zeta_zero(1).imag > 10
    def evaluate(poly, x):
        acc = arb(0)
        for coefficient in reversed(poly.coeffs()):
            acc = acc*x + arb(coefficient)
        return acc
    def xcoord(n):
        t = acb.zeta_zero(n).imag
        return -1/(2*(t*t+arb(1)/4))
    xs = [xcoord(n) for n in (11430, 11431, 11432)]
    assert xs[0] < xs[1] < xs[2]
    domain = arb(-1).union(arb(0))/200
    assert evaluate(denom, domain) > 0
    gaps = []
    for num in nums:
        derivative = num.derivative()*denom-num*denom.derivative()
        slope = evaluate(derivative, domain)
        assert slope > 0 or slope < 0
        readings = [evaluate(num, x)/evaluate(denom, x) for x in xs]
        left, right = [abs(readings[i]-readings[1])/2 for i in (0, 2)]
        assert right < left
        gaps.append(right)
    best = max(range(19), key=lambda k: float(gaps[k]))
    tied = [j for j in range(19) if nums[j] == nums[best]]
    assert all(gaps[j] < gaps[best] for j in range(19) if j not in tied)
    print("Single-port maximum attaining sites", tied, "critical absolute radius",
          gaps[best].str(21), flush=True)


if __name__ == "__main__":
    verify()
