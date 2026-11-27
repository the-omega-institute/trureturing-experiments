"""在两壳的四个精确等式上，检验第九壳牛顿与真空读数的局部独立性。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    ordered = lib("two-shell-ordered-bond-gates")
    box, _ = ordered["experiment"](root, fifth, columns, levels)
    active = ordered["ACTIVE"] + (3, 5)
    fixed = [arb(ordered["ANGLES"][i]) for i in (3, 5)]
    shells = []
    for number in range(2_306_888_100_092, 2_306_888_100_095):
        _, _, _, basis, _, weights, heat, field = lib("fifth-seventh-shortest-site-gates")["geometry"](
            root, fifth, columns, levels, number=number
        )
        shells.append((basis.conjugate().transpose(), weights, heat, field))
    gravity = lib("shortest-common-complex-gate")["gravity_system"]

    def matrix(point):
        gates = ordered["transport"](fifth, point, active=active)
        rows = []
        for left, weights, heat, field in shells:
            proj, *directions = [(left * gate) * columns for gate in gates]
            _, jac = gravity(proj, directions, levels, weights, heat, field, False)
            rows.extend(jac.tolist())
        return arb_mat(rows)

    center = [v.mid() for v in box]
    inverse = matrix(center + fixed).mid().inv()
    identity = arb_mat([[int(i == j) for j in range(6)] for i in range(6)])
    # 认证根在原盒中；将六个角再各向外放宽，以容纳以真根为中心的小盒。
    radius = arb("1e-6")
    angle_error = arb(0, radius.upper())
    stack = [[q + angle_error for q in box]]
    free = [q + angle_error for q in fixed]
    count = 0
    worst = arb(0)
    while stack:
        cell = stack.pop()
        defect = identity - inverse*matrix(cell + free)
        norm = max(sum(v.abs_upper() for v in row) for row in defect.tolist())
        if norm < arb(".9"):
            count += 1
            worst = max(worst, norm)
            continue
        axis = max(range(4), key=lambda i: cell[i].rad().upper())
        x = cell[axis]
        assert x.rad() > arb("1e-7")
        mid = x.mid()
        stack.extend([cell[:axis]+[edge.union(mid)]+cell[axis+1:]
                      for edge in (x.lower(), x.upper())])
    print("覆盖根盒的子盒数与六维雅可比的最大预条件偏差", count, worst.str(14), flush=True)
    gain = max(abs(inverse[i, 4]).upper() + abs(inverse[i, 5]).upper() for i in range(6))
    assert gain > 0
    output_radius = radius*(1-worst)/gain
    print("第九壳两个读数可独立指定的共同半径下界", output_radius.lower().str(15), flush=True)
    assert output_radius > arb("8.80e-9")
    return worst, output_radius


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))
