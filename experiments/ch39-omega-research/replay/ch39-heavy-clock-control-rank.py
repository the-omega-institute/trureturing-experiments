"""两只钟和径向亏额独立变化时，从目标塔的前两谱矩确定电导控制量的维数下界。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def native_mass_derivatives(shell, charge, height, deficit, radial, half_width):
    """在两角和亏额的参数盒内以 Hellmann–Feynman 公式计算十九个质量梯度。"""
    turn = acb(0, arb(0, half_width)).exp()
    phases = {d: z * turn for d, z in shell["geo"]["phases"](height).items()}
    native = shell["laplacian"](charge, height, phases=phases)
    dx, dy, dr = (acb_mat(19, 19) for _ in range(3))
    pole = shell["soft"]["pole"]
    tuned = shell["soft"]["star"]["orbit"]
    for i, j, distance, _ in shell["geo"]["edges"]:
        if pole in (i, j):
            if i in tuned or j in tuned:
                dr[i, i] -= 1
                dr[j, j] -= 1
                dr[i, j] += 1
                dr[j, i] += 1
        elif distance in (724, 804):
            derivative = dx if distance == 724 else dy
            derivative[i, j] = acb(0, charge) * native[i, j]
            derivative[j, i] = derivative[i, j].conjugate()
    native = native + dr * (radial - deficit)
    center = 15 - radial / 3
    dc = (arb(0), arb(0), -arb(1) / 3)
    values, vectors = native.eig(right=True, algorithm="rump")
    assert len(values) == 19 and all(v.imag.contains(0) and v.real > 0 for v in values)
    eigenvalues = sorted((v.real for v in values), key=float)
    assert all(a < b for a, b in zip(eigenvalues, eigenvalues[1:]))
    adjoint = vectors.conjugate().transpose()
    squared_norms = adjoint * vectors
    forms = [adjoint * derivative * vectors for derivative in (dx, dy, dr)]
    levels = []
    for k, value in enumerate(values):
        lam = value.real
        assert lam < center or lam > center
        sign = 1 if lam > center else -1
        norm = squared_norms[k, k].real
        assert norm > 0 and all(form[k, k].imag.contains(0) for form in forms)
        gradient = tuple(sign * ((2 * lam - center) * form[k, k].real / norm - lam * d)
                         for form, d in zip(forms, dc))
        levels.append((lam * abs(lam - center), gradient))
    levels.sort(key=lambda level: float(level[0]))
    assert all(a[0] < b[0] for a, b in zip(levels, levels[1:]))
    assert levels[2][0] < 10 < levels[3][0]
    return center, dc, levels


def moment_branches(center, dc, levels):
    """枚举轻质量的三个反解，返回每个分支的迹值及前两矩的三分量梯度。"""
    def root_and_gradient(level, scale, sign, side):
        mass, gradient = level
        discriminant = center**2 + 4 * sign * scale * mass
        assert discriminant > 0
        root = (center + side * discriminant.sqrt()) / 2
        assert root > 0
        derivative = tuple((root * d + sign * scale * g) / (2 * root - center)
                           for g, d in zip(gradient, dc))
        return root, derivative

    tolerance = arb("1e-8")
    assert 32 * levels[3][0] - tolerance > center**2 / 4
    assert center > 14 and (center**2 - 4 * (levels[2][0] + tolerance)).sqrt() > 13
    assert center < 2 * 13
    heavy = [root_and_gradient(level, 32, 1, 1) for level in levels[3:]]
    lights = [[root_and_gradient(level, 1, sign, side)
               for sign, side in ((1, 1), (-1, 1), (-1, -1))] for level in levels[:3]]
    branches = []
    for light in product(*lights):
        roots = heavy + list(light)
        trace = sum(z for z, _ in roots)
        first = tuple(sum(g[j] for _, g in roots) for j in range(3))
        second = tuple(sum(2 * z * g[j] for z, g in roots) for j in range(3))
        branches.append((trace, first, second))
    traces = sorted((t for t, _, _ in branches), key=float)
    assert all(b - a > 2 * 19 * tolerance / 13 for a, b in zip(traces, traces[1:]))
    return branches


def area(a, b):
    return a[0] * b[1] - a[1] * b[0]


def volume(a, b, c):
    return a[2] * area(b, c) + b[2] * area(c, a) + c[2] * area(a, b)


def verify():
    shell = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))
    ctx.prec = 320
    height = acb.zeta_zero(2_306_888_100_085).imag
    deficit = arb(23) / 25 * arb(2)**-20 / (height**2 + shell["soft"]["point"]["study"]["offset"])
    assert 0 < deficit < arb("1e-12")
    radial = arb(0).union(arb("1e-12"))
    half_width = arb("1e-9")
    branches = {}
    for q in (2, 4, 6):
        branches[q] = moment_branches(*native_mass_derivatives(shell, q, height, deficit, radial, half_width))
        floor = min(abs(area(first, second)).lower() for _, first, second in branches[q])
        assert floor > {2: 172, 4: 2090, 6: 6570}[q]
        # 任一线段的平均雅可比仍在同一逐项区间内。最小奇异值 ≥ |det|/Frobenius 范数。
        stretch = min((abs(area(a, b)) / sum(v*v for v in (*a[:2], *b[:2])).sqrt()).lower()
                      for _, a, b in branches[q])
        error = half_width * stretch / 58
        assert error > arb({2: "5.89e-12", 4: "1.95e-11", 6: "3.14e-11"}[q])
        print("荷、单塔两矩的角向行列式下界、单标量的逐层绝对误差下界",
              q, floor.str(12), error.str(12))
    for q, r in combinations(branches, 2):
        pairs = list(product(branches[q], branches[r]))
        angle_floor = min(abs(area(a, b)).lower() for (_, a, _), (_, b, _) in pairs)
        full_floor = min(abs(volume(a, a2, b)).lower() for (_, a, a2), (_, b, _) in pairs)
        assert angle_floor > arb({(2, 4): ".217", (2, 6): "2.15", (4, 6): ".267"}[q, r])
        assert full_floor > {(2, 4): 1408, (2, 6): 5937, (4, 6): 4543}[q, r]
        print("两荷、两迹的角向下界、前一荷两矩和后一荷迹的三维下界",
              q, r, angle_floor.str(12), full_floor.str(12))
    # 三个一阶迹单独成坐标的行列式较小，改在同心的小角盒上认证。
    branches = {q: moment_branches(*native_mass_derivatives(shell, q, height, deficit, radial,
                                                           arb("1e-13"))) for q in (2, 4, 6)}
    triple_floor = min(abs(volume(a, b, c)).lower() for (_, a, _), (_, b, _), (_, c, _)
                       in product(branches[2], branches[4], branches[6]))
    assert triple_floor > arb(".000239")
    print("三荷全部19683种分支的三迹体积下界", triple_floor.str(12))


if __name__ == "__main__":
    verify()
