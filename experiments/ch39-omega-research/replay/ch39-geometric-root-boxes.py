"""四个实际零点排除无钟子图给出的六根十二次作用及其独立扰动。

每部门选一个含三站轨道的 C3 不变子集，再选诱导图或原拉普拉斯主子块。
两个部门的六个字符最低本征值各可独立移动 1/200。对整个六维闭盒，
比较原十九站壳上的作用和 Σ_j Π_i (e_j-r_i)²，差额方向取目标减竞争者。
"""
from collections import Counter
from itertools import product
from math import prod
from pathlib import Path
from runpy import run_path
from typing import NamedTuple

import numpy as np
from flint import acb, acb_mat, arb, ctx


REPLAY = Path(__file__).resolve().parent
TARGETS = {-4: 2345666, 2: 964890}
RIVALS = ((-4, 111), (-4, 83), (2, 33), (2, 1))
SCALE = 10**9
HALF_WIDTH = SCALE // 200
ACTION_SCALE = SCALE**12


def lowest_on_grid(matrix, adjoint):
    """以一格为误差，认证 Hermitian 块的最低本征值，可含重根。"""
    approx = np.array([[complex(z) for z in row] for row in matrix.tolist()])
    levels, vectors = np.linalg.eigh(approx)
    center = round(float(levels[0]) * SCALE)
    lower = arb(center - 1) / SCALE

    # Sylvester 判据给出严格下界。浮点向量只用于提出 Rayleigh 商的候选。
    for size in range(1, matrix.nrows() + 1):
        minor = acb_mat([
            [matrix[i, j] - (lower if i == j else 0) for j in range(size)]
            for i in range(size)
        ])
        assert minor.det().real > 0
    vector = acb_mat([[acb(complex(z))] for z in vectors[:, 0]])
    norm_squared = (adjoint(vector) * vector)[0, 0].real
    numerator = (adjoint(vector) * matrix * vector)[0, 0].real
    assert norm_squared > 0 and numerator / norm_squared < arb(center + 1) / SCALE
    return center


def reference_triples(shell, pump, charge):
    """按子集和两种边界约定，给出 252 组三字符参考根的整数包围中心。"""
    adjoint = pump["adjoint"]
    frames = list(pump["rotation_bases"](shell))
    base, links = pump["clock_matrices"](shell, charge)
    link = links[0] + links[1]
    full = base + link + adjoint(link)
    edges = [(i, j) for i, j, _, _ in shell["geo"]["edges"]]

    # 每个平凡字符基列恰好支撑在一个轨道上，非零分量的模为 1 或 1/√3。
    trivial = frames[0]
    orbits = [[i for i in range(trivial.nrows()) if abs(trivial[i, j]) > arb(".1")]
              for j in range(trivial.ncols())]
    assert sorted(map(len, orbits)) == [1, 3, 3, 3, 3, 3, 3]
    triples = []
    for mask in range(1, 1 << len(orbits)):
        sites = sorted(i for j, orbit in enumerate(orbits) if mask & (1 << j) for i in orbit)
        if len(sites) == 1:
            continue
        degrees = dict.fromkeys(sites, 0)
        for i, j in edges:
            if i in degrees and j in degrees:
                degrees[i] += 1
                degrees[j] += 1
        principal = acb_mat([[full[i, j] for j in sites] for i in sites])
        induced = acb_mat(principal)
        for i, site in enumerate(sites):
            induced[i, i] = degrees[site]

        restricted_frames = []
        for frame in frames:
            columns = [k for k in range(frame.ncols()) if any(abs(frame[i, k]) > arb(".1") for i in sites)]
            restricted_frames.append(acb_mat([[frame[i, k] for k in columns] for i in sites]))
        for matrix in (induced, principal):
            triples.append([lowest_on_grid(adjoint(basis) * matrix * basis, adjoint)
                            for basis in restricted_frames])
    assert len(triples) == 252
    return triples


def actual_levels(light, charge, index):
    centers = []
    for energy in light(charge, index):
        center = round(float(energy * SCALE))
        assert abs(energy * SCALE - center) < 1
        centers.append(center)
    return centers


class FactorBounds(NamedTuple):
    lower: int
    upper: int
    derivative: int
    bernstein: np.ndarray


def triple_bounds(roots, energy):
    """包围三个平方因子的乘积。所有量使用整数格单位。"""
    # 一格覆盖参考本征值的取整误差，再加 1/200 的独立扰动。
    radius = HALF_WIDTH + 1
    distances = [(energy - root - radius, energy - root + radius) for root in roots]
    # 再把实际能量放宽一格。最大绝对距离也用于能量导数的乘积法则。
    maxima = [max(abs(a), abs(b)) + 1 for a, b in distances]
    lower = prod(max(0, a - 1, -b - 1) for a, b in distances)**2
    upper = prod(maxima)**2
    derivative = sum(2 * d * prod(v**2 for k, v in enumerate(maxima) if j != k)
                     for j, d in enumerate(maxima))

    # 当 e-r 从 a 走到 b，平方的二次 Bernstein 系数为 (a²,ab,b²)。
    # 三根独立变化，27 个张量积系数均用 Python 任意精度整数保存。
    factors = [(a*a, a*b, b*b) for a, b in distances]
    coefficients = np.array([prod(terms) for terms in product(*factors)], dtype=object)
    return FactorBounds(lower, upper, derivative, coefficients)


def interval_margin(up, down, target, rival):
    return (sum(up[e].lower * down[e].lower for e in target)
            - sum(up[e].upper * down[e].upper for e in rival))


def bernstein_margin(up, down, target, rival):
    # 六根的张量积 Bernstein 基函数非负且和为一，最小系数给出全盒下界。
    coefficients = np.zeros((27, 27), dtype=object)
    rounding_error = 0
    for sign, levels in ((1, target), (-1, rival)):
        for energy in levels:
            left, right = up[energy], down[energy]
            coefficients += sign * np.outer(left.bernstein, right.bernstein)
            # 以整数能量代替真实能量时相差至多一格，乘积导数包住整段。
            rounding_error += left.derivative * right.upper + left.upper * right.derivative
    return min(coefficients.flat) - rounding_error


def verify():
    shell = run_path(str(REPLAY / "ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(REPLAY / "ch39-native-light-projector-pump.py"))
    selector = run_path(str(REPLAY / "ch39-light-spectrum-zero-selector.py"))
    light = selector["light_spectra"](REPLAY)
    ctx.prec = 256
    targets = {q: actual_levels(light, q, n) for q, n in TARGETS.items()}
    comparisons = [(targets[q], actual_levels(light, q, n)) for q, n in RIVALS]
    energies = {e for pair in comparisons for levels in pair for e in levels}
    tables = []
    for charge in (-4, 2):
        tables.append([{e: triple_bounds(roots, e) for e in energies}
                       for roots in reference_triples(shell, pump, charge)])

    counts = Counter()
    refinements = 0
    minimum = None
    for up, down in product(*tables):
        margins = [interval_margin(up, down, target, rival) for target, rival in comparisons]
        order = sorted(range(len(RIVALS)), key=margins.__getitem__, reverse=True)
        witness = order[0]
        bound = margins[witness]
        # 乘积区间给出百单位以下的余量时，按同一顺序检查更紧的系数界。
        if bound <= 100 * ACTION_SCALE:
            refinements += 1
            for witness in order:
                bound = bernstein_margin(up, down, *comparisons[witness])
                if bound > 100 * ACTION_SCALE:
                    break
        assert bound > 112 * ACTION_SCALE
        counts[RIVALS[witness]] += 1
        minimum = bound if minimum is None else min(minimum, bound)

    total = sum(counts.values())
    assert total == 252**2
    print("六根独立扰动半宽 1/200，配方闭盒数", total)
    print("乘积区间与 Bernstein 系数分别认证", total - refinements, refinements)
    print("荷、竞争编号及其认证的闭盒数", [(q, n, counts[q, n]) for q, n in RIVALS])
    # 十进制向下取整输出，保留可直接引用的严格余量。
    places = 10**10
    whole, decimal = divmod(minimum * places // ACTION_SCALE, places)
    print(f"目标减竞争者的统一下界 {whole}.{decimal:010d}")


if __name__ == "__main__":
    verify()
