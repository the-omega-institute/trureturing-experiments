"""固定第 11661 枚零点的六面味坐标时，连续双钟相对原全局真空的最小作用代价。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def circle_minimum(model, modes, gap_angle, center):
    # 固定相角差后 y=t、x=t+d。第三个谐波与 t 无关。
    def action(t):
        return sum(model["energy"](q, coefficients, t+gap_angle, t)
                   for q, coefficients in modes)

    def derivatives(t):
        slope = arb(0)
        curvature = arb(0)
        for q, (a, b, _) in modes:
            value = (a * acb(0, q*gap_angle).exp() + b) * acb(0, q*t).exp()
            slope += q * value.imag
            curvature += q*q * value.real
        return slope, curvature

    radius = arb(".002")
    trial = action(center)
    pending = [(1, 0)]
    while pending:
        n, i = pending.pop()
        t = arb.pi() * (arb(2*i+1, 1)/n-1)
        if abs(t-center) < radius or action(t) > trial:
            continue
        slope, _ = derivatives(t)
        if not slope.contains(0):
            continue
        assert n < 2**19, (n, i, t)
        pending.extend(((2*n, 2*i), (2*n, 2*i+1)))

    _, curvature = derivatives(center + arb(0, radius))
    assert curvature > 100
    for _ in range(4):
        slope, second = derivatives(center)
        center = (center - slope/second).mid()
    assert abs(center - (arb(".11101666") if gap_angle < 0 else arb("-.11101666"))) < radius/2
    slope, _ = derivatives(center)
    enclosure = center + arb(0, slope.abs_upper()/100)
    return action(enclosure), enclosure


def three_window_gap(root, local, model, modes, angles, ground):
    # 原曲率必要谱带包含对三条软辐边的扰动余量。
    window = run_path(str(root / "ch39-curvature-window-density-cap.py"))
    lower, upper = window["eigen_band"]()
    below = local["spectral_cut"](root, lower.lower())
    above = local["spectral_cut"](root, upper.upper())
    for angle in angles:
        pending = [(1, 0)]
        while pending:
            n, i = pending.pop()
            y = arb.pi() * (arb(2*i+1, 1)/n - 1)
            x = y + angle
            value = sum(model["energy"](q, c, x, y) for q, c in modes)
            if value - ground > arb("1.3266175"):
                continue
            if below(x, y, angle)[0] >= 16 or above(x, y, angle)[1] < 16:
                continue
            assert n < 2**35, (angle, y, value-ground)
            pending.extend(((2*n, 2*i), (2*n, 2*i+1)))

    # 在曲率窗边界内侧给出同一高度范围通用的可行上界。
    y = arb("-.13699465")
    x = y + angles[1]
    physical = run_path(str(root / "ch39-first-three-mode-observables.py"))
    local["full_windows"](physical, x, y)
    trial_gap = sum(model["energy"](q, c, x, y) for q, c in modes) - ground
    assert trial_gap < arb("1.32661933")
    print("保持三窗和参考味的能隙下界、可行见证能隙", "1.3266175", trial_gap.str(17))


def verify():
    root = Path(__file__).resolve().parent
    model = run_path(str(root / "ch39-mincut-global-winding.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    flavor = run_path(str(root / "ch39-chiral-neutral-flavor-phase.py"))
    ctx.prec = 224
    tri = model["terms"](shell)
    assert sum(map(len, tri.values())) == 529 and shell["soft"]["point"]["study"]["offset"] > 0
    delta = arb(0).union(arb(23)/25*arb(2)**-20/10**24)
    modes = [(q, [w*c/113 for c in model["coeff"](tri, q, delta)])
             for q, w in ((1, 100), (27, 13))]
    center = model["global_cover"](tri)
    radius = arb(0, "2.1e-17")
    vacuum = [sum(model["energy"](q, c, sign*center[0]+radius, sign*center[1]+radius)
                  for q, c in modes) for sign in (-1, 1)]
    ground = vacuum[0].union(vacuum[1])
    local = run_path(str(root / "ch39-mincut-reference-ckm-action.py"))
    angles = local["reference_angles"](flavor)
    for gap_angle in angles:
        center_y = arb(".11101666") if gap_angle < 0 else arb("-.11101666")
        best, y = circle_minimum(model, modes, gap_angle, center_y)
        gap = best - ground
        assert arb("1.25008749932364") < gap < arb("1.25008749932367")
        print("参考相角差、约束极小点第二钟角、归一作用能隙", gap_angle.str(17),
              y.str(20), gap.str(17))
    three_window_gap(root, local, model, modes, angles, ground)


if __name__ == "__main__":
    verify()
