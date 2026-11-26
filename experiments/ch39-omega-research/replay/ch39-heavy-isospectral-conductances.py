"""固定十九层的目标质量，在逐条边的116个独立坐标上同时构造局部等谱族。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def max_incidence(edges):
    return max(sum(i in edge for edge in edges) for i in range(19))


def verify():
    root = Path(__file__).resolve().parent
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    towers = run_path(str(root / "ch39-three-even-positive-towers.py"))
    continuation = run_path(str(root / "ch39-three-tower-phase-continuation.py"))
    ctx.prec = 320
    height = acb.zeta_zero(towers["NUMBER"]).imag
    delta = arb(23)/25 * arb(2)**-20 / (height**2 + shell["soft"]["point"]["study"]["offset"])
    center = 15 - delta/3
    groups = towers["orbits"](shell["geo"])
    free_width, solved_width = arb("1e-13"), arb("1e-9")
    for q in (2, 4, 6):
        native = shell["laplacian"](q, height)
        num, variable, perm = towers["DATA"][q]
        w = [arb(n)/10**15 for n in num]
        phases = {(i, j): -native[i, j] / abs(native[i, j]) for group in groups for i, j in group}
        representatives = [groups[g][0] for g in variable]
        remainder = sorted(set(phases) - set(representatives))
        assert (len(representatives), len(remainder)) == (19, 116)
        raw = native.eig(multiple=True)
        assert len(raw) == 19 and all(v.imag.contains(0) and v.real > 0 for v in raw)
        folded = sorted((v.real * abs(v.real - center) for v in raw), key=float)
        assert all(a < b for a, b in zip(folded, folded[1:]))
        target, _, _ = continuation["targets"](folded, center, perm)
        assert all(a < b for a, b in zip(target, target[1:]))
        matrix = acb_mat(19, 19)
        for g, group in enumerate(groups):
            for i, j in group:
                matrix[i, i] += w[g]
                matrix[j, j] += w[g]
                matrix[i, j] = -w[g] * phases[i, j]
                matrix[j, i] = matrix[i, j].conjugate()
        eig, vec = matrix.eig(right=True, algorithm="rump")
        order = sorted(range(19), key=lambda i: float(eig[i].real))
        levels = [eig[k].real for k in order]
        gap = min((b-a).lower() for a, b in zip(levels, levels[1:]))
        derivative = arb_mat(19, 19)
        for a, k in enumerate(order):
            norm = sum(abs(vec[i, k])**2 for i in range(19))
            for b, (i, j) in enumerate(representatives):
                derivative[a, b] = abs(vec[i, k]-phases[i, j]*vec[j, k])**2/norm
        pre = derivative.inv()
        pre_norm = continuation["rownorm"](pre)
        correction = continuation["rownorm"](pre * arb_mat([[a-b] for a, b in zip(levels, target)]))
        free_perturbation = 2 * max_incidence(remainder) * free_width
        change = free_perturbation + 2 * max_incidence(representatives) * solved_width
        assert gap > 2 * change
        # 每条边的拉普拉斯导数范数为2，单根投影迹范数的变化至多2b/(gap-2b)。
        contraction = 19 * pre_norm * 4 * change / (gap-2*change)
        residual = correction + pre_norm * free_perturbation
        assert contraction < 1 and residual + contraction * solved_width < solved_width
        assert min(w) > solved_width + free_width
        print("荷、逐边逆雅可比范数、预乘残差上界、整盒压缩率",
              q, pre_norm.str(12), residual.str(12), contraction.str(12))


if __name__ == "__main__":
    verify()
