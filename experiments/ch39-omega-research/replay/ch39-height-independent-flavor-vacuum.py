"""同一组正绕数权在全部正零点高度将同一精确参考味相位选为全局真空。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def fixed_weights(lib, tri, x, y):
    # 八个有理数固定一次，余下五权仅由不含高度的五个线性方程决定。
    fixed = {1: ".0114566", 3: ".0132862", 28: ".0766797", 33: ".0506797",
             40: ".0212736", 69: ".0891141", 112: ".0228542", 125: ".0018878"}
    return solve_weights(lib, tri, x, y, fixed, [7, 35, 42, 83, 105])


def minimum_order_weights(lib, tri, x, y):
    fixed = {1: ".0198450", 2: ".0067192", 3: ".0238946", 4: ".0115864", 28: ".0107808",
             65: ".0010449", 67: ".0197454", 83: ".0113660", 104: ".0407273"}
    return solve_weights(lib, tri, x, y, fixed, [30, 35, 38, 40, 105])


def gradient_columns(lib, tri, q, x, y):
    bare = lib["coeff"](tri, q, arb(0))
    soft = [sum(z**q for z, count in tri[k] if count == 1)
            for k in ((1, 0), (0, 1), (1, -1))]
    return [*lib["derivatives"]([(q, bare)], x, y)[:2],
            *lib["derivatives"]([(q, soft)], x, y)[:2]]


def solve_weights(lib, tri, x, y, fixed, unknown):
    fixed = {q: arb(w) for q, w in fixed.items()}
    columns = {}
    for q in [*fixed, *unknown]:
        columns[q] = [arb(1), *gradient_columns(lib, tri, q, x, y)]
    target = [arb(int(k == 0))-sum(w*columns[q][k] for q, w in fixed.items())
              for k in range(5)]
    matrix = arb_mat([columns[q] for q in unknown]).transpose()
    assert not matrix.det().contains(0)
    solution = matrix.solve(arb_mat([[v] for v in target])).entries()
    weights = {**fixed, **dict(zip(unknown, solution, strict=True))}
    assert min(weights.values()) > arb(".001")
    return weights


def exclude_lower_orders(lib, tri, x, y):
    # 全部实际高度全局取小，经连续性推出单位辐边极限同样全局取小。
    # 三组比较点只在单位辐边上评能。两个不同高度则给出四个精确梯度等式。
    comparison = [("-.0073015873", "-.1387301587", ".10128652"),
                  (".0073015873", "-.1387301587", "277.53755061"),
                  (".0133333333", "-.0488888889", "3.12878691")]
    directions = [arb(v) for v in ("1.14312953", "-.07014091", ".32261739", ".16782987")]
    margins = []
    for q in range(1, 105):
        bare = lib["coeff"](tri, q, arb(0))
        original = lib["energy"](q, bare, x, y)
        margin = sum(arb(w)*(original-lib["energy"](q, bare, arb(a), arb(b)))
                     for a, b, w in comparison)
        margin += sum(u*v for u, v in zip(directions, gradient_columns(lib, tri, q, x, y), strict=True))
        margins.append(margin)
    assert min(margins) > arb(".99")
    print("同一相位全高度取小：前104阶逐阶严格分离的最小余量", min(margins).str(13))


def required_high_weights(lib, tri, x, y):
    # 比较点坐标是 1/56 的整数倍。权重与乘子写为给定十进有理数。
    certificates = {
        104: (
            [(-1, -9, '0.02741545404'), (-1, -8, '0.02415973638'), (-2, -7, '0.00909936582'), (1, -3, '0.01403918796'), (2, -3, '0.01774241032'), (1, -1, '0.00028141381')],
            ['-0.00273319126', '-0.00060654571', '-0.00152459297', '0.00010813929'],
            '.006589'),
        105: (
            [(-1, -8, '0.09842425062'), (1, -3, '0.07891318210'), (1, -2, '0.14083802372'), (1, -1, '0.03966768135')],
            ['-0.01198210623', '-0.00116425477', '-0.00624742343', '-0.00197606174'],
            '.068674'),
    }
    for order, (points, multipliers, bound) in certificates.items():
        directions = [arb(u) for u in multipliers]
        margins = []
        for q in range(1, 106):
            bare = lib["coeff"](tri, q, arb(0))
            trial = lib["energy"](q, bare, x, y)
            margin = arb(int(q == order))
            for a, b, w in points:
                rival = lib["energy"](q, bare, arb(a)/56, arb(b)/56)
                margin += arb(w)*(trial-rival)
            margin -= sum(t*v for t, v in zip(directions, gradient_columns(lib, tri, q, x, y), strict=True))
            margins.append(margin)
        assert min(margins) > arb(bound)
        print("刚达到105圈时不可缺少的阶数和归一权下界", order, min(margins).str(12))


def verify():
    root = Path(__file__).resolve().parent
    lib = run_path(str(root / "ch39-mincut-global-winding.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    global_model = run_path(str(root / "ch39-reference-flavor-global-winding.py"))
    local = run_path(str(root / "ch39-mincut-reference-ckm-action.py"))
    flavor = run_path(str(root / "ch39-chiral-neutral-flavor-phase.py"))
    ctx.prec = 224
    tri = lib["terms"](shell)
    global_model["conjugation_symmetry"](root)
    cap = arb(23)/25*arb(2)**-20/14**2
    assert acb.zeta_zero(1).imag > 14 and shell["soft"]["point"]["study"]["offset"] > 0
    x = arb("-.01046")
    y = x-local["reference_angles"](flavor)[0]
    exclude_lower_orders(lib, tri, x, y)
    required_high_weights(lib, tri, x, y)
    delta = arb(0).union(cap)
    for weight_rule in (fixed_weights, minimum_order_weights):
        weights = weight_rule(lib, tri, x, y)
        modes = [(q, [w*z for z in lib["coeff"](tri, q, delta)]) for q, w in weights.items()]
        drift = global_model["global_cover"](lib, tri, weights, modes, x, y)
        assert drift < arb("-.05")
    global_model["physical_all_heights"](root, x, y, cap)


if __name__ == "__main__":
    verify()
