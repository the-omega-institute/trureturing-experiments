"""用三个荷的质量平方和检验半周拼接能否由同一谱权产生。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx

shell = run_path(str(Path(__file__).with_name("ch39-monotone-heavy-mass-ratio.py")))
physical = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))
charges = (1, 3, 2)
frequencies = range(-4, 5)
grid_size = 256


def moment_coefficients(charge):
    # 在 z^q、w^q 平面取样。L 的矩阵元一次，平方质量迹至多四次。
    ninth = [acb(0, 2 * arb.pi() * i / 9).exp() for i in range(9)]
    cut = shell["cut_matrix"](charge)
    edges = [(i, j, d, -u**charge) for i, j, d, u in shell["SHELL"]["edges"]
             if d in (724, 804)]
    values = acb_mat(9, 9)
    for a, z in enumerate(ninth):
        for b, w in enumerate(ninth):
            matrix = acb_mat(cut)
            for i, j, d, value in edges:
                matrix[i, j] = value * (z if d == 724 else w)
                matrix[j, i] = matrix[i, j].conjugate()
            product = matrix * matrix - 15 * matrix
            # 在单位相位上矩阵厄米，质量平方和等于该矩阵的平方迹。
            values[a, b] = (product * product).trace().real
    analysis = acb_mat([[ninth[-k*i % 9] for i in range(9)] for k in frequencies])
    return analysis * values * analysis.transpose() / 81


def moment_at(coefficients, z, w):
    row = acb_mat([[z**k for k in frequencies]])
    column = acb_mat([[w**k] for k in frequencies])
    return (row * coefficients * column)[0, 0].real


def grid_enclosure(q, coefficients):
    roots = [acb(0, 2 * arb.pi() * i / grid_size).exp() for i in range(grid_size)]
    evaluation = acb_mat([[roots[i*k*q % grid_size] for k in frequencies]
                          for i in range(grid_size)])
    values = evaluation * coefficients * evaluation.transpose()
    first = arb(0)
    second = arb(0)
    for i, a in enumerate(frequencies):
        for j, b in enumerate(frequencies):
            size = abs(coefficients[i, j])
            first += q * (abs(a) + abs(b)) * size
            second += q**2 * (a*a + b*b) * size
    # 双线性插值余项，两个方向各有 h²/8 的二阶导数界。
    return values, second * (2 * arb.pi() / grid_size)**2 / 8, first


def verify():
    ctx.prec = 160
    data = {q: moment_coefficients(q) for q in charges}
    grids = {q: grid_enclosure(q, data[q]) for q in charges}
    # 真实辐边的谱移至多 4δ，中性中心移动 δ/3，全部正零点高度超过 14。
    delta = arb(23) / 25 * arb(2)**-20 / 14**2
    soft_error = 19 * 1600 * (4 * 80 + 12) * delta
    assert soft_error < arb(".05")
    window = arb(2)**-32
    z = acb(0, arb.pi()/64).exp()
    w = z**2
    widths = {}
    for q in charges:
        values, error, _ = grids[q]
        samples = [v.real for v in values.entries()]
        widths[q] = max(v.upper() for v in samples) - min(v.lower() for v in samples)
        widths[q] += 2 * (error + soft_error)

    witnesses = {}
    for n in (10, 26, 33, 36, 127, 200, 315, 316, 946):
        height = acb.zeta_zero(n).imag
        assert height > 14
        phases = shell["SHELL"]["phases"](height)
        center = physical["laplacian"](0, height, phases=phases).trace().real / 18
        witnesses[n] = {}
        for q in charges:
            laplacian = physical["laplacian"](q, height, phases=phases)
            product = laplacian * laplacian - center * laplacian
            value = (product * product).trace().real
            assert abs(value - moment_at(data[q], phases[724]**q, phases[804]**q)) < soft_error
            witnesses[n][q] = value

    for shifted in (1, 3):
        scale = {q: arb(8)/11 if q == shifted else arb(-1) for q in charges}
        target = {}
        for q in charges:
            phase = -w if q == shifted else w
            target[q] = moment_at(data[q], z**q, phase**q)
            target[q] += arb(0, soft_error + window * grids[q][2])
        for pair in ((1, 3), (1, 2), (3, 2)):
            if pair == (1, 3):
                indices = (10, 26, 33)
            elif pair[0] != shifted:
                indices = (315, 316, 946)
            else:
                indices = (33, 36 if shifted == 1 else 127, 200)
            vertices = arb_mat([[witnesses[n][q] for n in indices] for q in pair] + [[1, 1, 1]])
            weights = vertices.solve(arb_mat([[target[q]] for q in pair] + [[1]]))
            assert all(v > arb(".015") for v in weights.entries())
            print("半周荷、配平的两种荷、零点编号与正权", shifted, pair, indices,
                  *[v.str(8) for v in weights.entries()])

        table = sum((grids[q][0] * scale[q] for q in charges),
                    acb_mat(grid_size, grid_size))
        maximum = max(v.real.upper() for v in table.entries())
        maximum += sum(abs(scale[q]) * (grids[q][1] + soft_error) for q in charges)
        margin = sum(scale[q] * target[q] for q in charges) - maximum
        # 若允许分荷选权，以荷二的分布为共同基准，极差乘总变差控制期望差。
        tv_bound = margin / sum(abs(scale[q]) * widths[q] for q in (1, 3))
        assert margin > (1022 if shifted == 1 else 1155)
        assert tv_bound > arb(".06459" if shifted == 1 else ".07254")
        print("转半周的荷、线性分离余量下界、分荷总变差下界", shifted,
              margin.lower().str(14), tv_bound.lower().str(14))


if __name__ == "__main__":
    verify()
