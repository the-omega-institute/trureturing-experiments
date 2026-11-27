"""用两个目的层的 Gram 矩阵核验两份自由补空间见证到三次单位根的距离。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def two_by_two_minimum(matrix):
    a, d = matrix[0, 0].real, matrix[1, 1].real
    gap = ((a - d)**2 + 4 * abs(matrix[0, 1])**2).sqrt()
    return (a + d - gap) / 2


def verify(root):
    stabilizers = run_path(str(root / "ch39-up-gravity-block-stabilizers.py"))
    spectral = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 384
    target = spectral["spectrum"](shell, -4, 2_306_888_100_087)[1]
    for sites, transform, limit in zip(stabilizers["BLOCKS"], stabilizers["build_inputs"](root),
                                        (".029", ".008"), strict=True):
        free = [j for j in range(19) if j not in sites]
        probes = acb_mat([[target[j, k] for k in (15, 17)] for j in free])
        block = acb_mat([[transform[i, j] for j in free] for i in free])
        gram = adj(probes) * probes
        cross = adj(probes) * block * probes
        # 首个三分之一圆周的网格及其全部旋转，单元中心到边缘至多 pi/(3n)。
        n = 1200
        correction = 2 * arb.pi() / (9 * n)
        margin = arb(limit)
        worst = []
        for i in range(n):
            minima = []
            for k in range(3):
                phase = acb(0, arb.pi() * (2*i+1+2*k*n) / (3*n)).exp()
                shifted = cross * phase
                matrix = (gram + shifted + adj(shifted)) * (arb(1)/3)
                minima.append(two_by_two_minimum(matrix).upper())
            worst.append(min(minima))
        upper = max(worst) + correction
        assert upper < -margin
        print("局部站、全部整体相位下最差投影的最小本征值上界、补空间距离下界",
              sites, upper.str(13), (3*margin/2).str(9), flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
