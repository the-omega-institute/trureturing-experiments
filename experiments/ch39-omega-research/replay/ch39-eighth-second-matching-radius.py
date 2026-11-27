"""第八壳上，按弦长共用站位门时，求恢复两项引力读数所需的最小距离。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    sites = lib("three-ray-chiral-boundary")["geo"]["sites"]
    pairs = {}
    for i, j in combinations(range(len(sites)), 2):
        length = sum((a-b)**2 for a, b in zip(sites[i], sites[j], strict=True))
        pairs.setdefault(length, []).append((i, j))
    assert sorted(pairs)[:2] == [12, 32]
    assert pairs[12] == [(6, 11), (7, 12), (8, 18)]
    assert pairs[32] == [(0, 2), (1, 3), (4, 5)]
    assert len({v for d in (12, 32) for e in pairs[d] for v in e}) == 12

    geometry = lib("fifth-seventh-shortest-site-gates")["geometry"]
    base, _, _, basis, station, weights, heat, field = geometry(
        root, fifth, columns, levels, number=2_306_888_100_093
    )
    # P_+ 和 P_- 的两条重列在第八壳本征标架的坐标。
    def projections(edges, sign):
        return acb_mat([[sum((basis[a, i] + sign*basis[b, i]).conjugate()
                             * (station[a, j] + sign*station[b, j]) / 2 for a, b in edges)
                         for j in range(2)] for i in range(19)])

    plus = [projections(pairs[d], 1) for d in (12, 32)]
    minus = [projections(pairs[d], -1) for d in (12, 32)]
    fixed = base - plus[0] - plus[1] - 2*(minus[0] + minus[1])

    def project(point):
        matrix = fixed
        directions = []
        for q, block in zip(point, plus, strict=True):
            phase = -acb(1-q**2, 2*q)/(1+q**2)
            matrix = matrix + phase*block
            directions.append(acb(0, 2)*phase*block/(1+q**2))
        return matrix, directions

    center = list(map(arb, (".38590088996362", "-.56884184746218")))
    certify = lib("shortest-common-complex-gate")["certify"]
    box, ratio, determinant = certify(project, center, arb(".0002"), levels, weights, heat, field)
    assert ratio < arb(".647") and not determinant.contains(0)
    return box


def verify(root):
    joint = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    return experiment(root, *joint["verify"](root, return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
