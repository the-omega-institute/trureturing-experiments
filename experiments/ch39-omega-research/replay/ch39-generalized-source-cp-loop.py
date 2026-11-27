"""用可测质量与 CKM 回路给出两部门广义本征基内的源回路相位。"""
from itertools import combinations
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb


def gap_ratio(spectrum):
    return prod(b-a for a, b in combinations(spectrum, 2)) / prod(spectrum)


def invariant(eigenvalues, readings):
    u, d = [[m**2 for m in sector] for sector in (readings[:3], readings[3:6])]
    t12, t23, t13, delta = readings[6:]
    jarlskog = t12.sin()*t12.cos()*t23.sin()*t23.cos()*t13.sin()*t13.cos()**2*delta.sin()
    return jarlskog * gap_ratio(u) * gap_ratio(d) / gap_ratio(eigenvalues)


def verify():
    def load(name):
        return run_path(str(Path(__file__).with_name(f"ch39-{name}.py")))

    model, a, b = load("two-zero-orbit-condition")["calibration"]()
    adj = model["adj"]
    relative = a * b.inv()
    eigenvalues, _ = model["eig"](adj(relative)*relative)
    readings = [x + arb(0, w.upper()) for x, w in zip(model["ref"], model["tol"], strict=True)]
    eta = invariant(eigenvalues, readings)
    assert 7400 < eta < 16300
    weights, tail = load("real-quadratic-flavor-window")["component_bounds"](model, relative, eigenvalues)
    magnitudes = []
    for a1, a2 in combinations(weights, 2):
        main = a1 * a2
        # 在最重下型行调实的标架中，每条边在以 main 为中心、tail 为半径的圆盘内。
        assert tail < main / 100
        magnitudes.append((main - tail).lower().union((main + tail).upper()))
    d, s, b_mass = readings[3:6]
    sine = eta * (d*s/b_mass**2)**2 / prod(magnitudes)
    # 每个圆盘的辐角绝对值 < asin(1/100)，三相位之和在反正弦主支。
    phase = sine.asin()
    assert arb(".00050") < phase < arb(".0057")
    print("广义源回路归一虚部、回路辐角下界与上界", eta.str(5), phase.lower().str(8), phase.upper().str(8))

    # 在实际零点上从协方差和十项读数两侧分别计算，并核对行列式恒等式的取向。
    gram = load("quadratic-gram-source")
    a, b = [gram["at_zero"](model, q, n) for q, n in ((-4, 2345666), (2, 964890))]
    relative = a*b.inv()
    eigenvalues, frame = model["eig"](adj(relative)*relative)
    k = model["K"]
    z = adj(frame)*b*k
    source = z*adj(z)
    loop = source[0, 1]*source[1, 2]*source[2, 0]
    direct = loop.imag / abs(z.det())**2
    assert 11519 < direct < 11521
    (mu, u), (md, d) = [model["eig"](adj(y)*y) for y in (a*k*arb("7.06539940386823251"), b*k)]
    obs = [v.sqrt() for v in mu+md] + model["mixing_angles"](u, d)
    from_observables = invariant(eigenvalues, obs)
    assert direct.overlaps(from_observables) and 11519 < from_observables < 11521
    assert arb(".001469") < loop.arg() < arb(".001470")
    print("原见证的两侧数值和源回路辐角", direct.str(10), from_observables.str(10), loop.arg().str(9))


if __name__ == "__main__":
    verify()
