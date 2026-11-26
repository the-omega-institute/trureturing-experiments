"""原壳正负荷双态的站点投影作用、实结构下的全角域界和实际零点。"""
from collections import Counter
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
cut = run_path(str(root / "ch39-monotone-heavy-mass-ratio.py"))
chiral = run_path(str(root / "ch39-chiral-common-gap-crossing.py"))
geometry = shell["geo"]
charges = (1, 2, 3, 4, 6)


def graph_and_triangle():
    sites = geometry["sites"]
    edges = {(i, j): (d, u) for i, j, d, u in geometry["edges"]}
    assert len(sites) == 19 and len(edges) == 135
    degrees = Counter(v for i, j in edges for v in (i, j))
    assert min(degrees.values()) == 10
    assert all(k*(11-k) >= 18 for k in range(2, 10))
    assert sum(d in (724, 804) for d, _ in edges.values()) == 12
    pole = shell["soft"]["pole"]
    tuned = shell["soft"]["star"]["orbit"]
    assert len(tuned) == 3 and degrees[pole] == 18 and all(degrees[v] == 10 for v in tuned)
    assert all((min(v, pole), max(v, pole)) in edges for v in tuned)

    assert sites[14] == (17, 17, 17) and sites[0] == (1, 5, 29) and sites[16] == (23, 7, 17)
    assert [edges[pair][0] for pair in ((0, 14), (0, 16), (14, 16))] == [544, 632, 136]
    z = edges[0, 14][1].conjugate()*edges[0, 16][1]*edges[14, 16][1].conjugate()
    exact = acb(2438004, -7208*arb(867).sqrt())
    assert (z-exact/abs(exact)).contains(0)
    assert arb(".08683") < abs(z.arg()) < arb(".08684")
    assert all(not (z**(2*q)).imag.contains(0) for q in charges)
    return z


def magnetic_matrix(q, clocks, *, omit_clock_edges=False):
    # 所有边先取单位权。证明中再统一乘最小实权平方以得到向下的算子界。
    out = acb_mat(19, 19)
    constant_trial = arb(0)
    for i, j, distance, transport in geometry["edges"]:
        if omit_clock_edges and distance in (724, 804):
            continue
        phase = (transport*clocks.get(distance, 1))**q
        out[i, i] += 1
        out[j, j] += 1
        out[i, j] = -phase**2
        out[j, i] = out[i, j].conjugate()
        constant_trial += 2*phase.imag**2
    return out, constant_trial


def inertia_lower_bound(matrix, min_weight_squared, energy):
    # 对 |z_v|=1 的十九维向量，E≥(19/2) min(w²) λ_min。
    level = 2*energy/(19*min_weight_squared)
    assert cut["count_above"](matrix, level) == 19


def verify():
    ctx.prec = 192
    triangle = graph_and_triangle()
    first = acb.zeta_zero(1).imag
    assert first > 14 and shell["soft"]["point"]["study"]["offset"] > 0
    delta_bound = arb(23)/25*arb(2)**-20/14**2
    edge_floor = (1-delta_bound)**2
    assert 18*edge_floor > 10 and 9+edge_floor > arb("9.99999999")
    floors = ("0.32331", "1.28531", "2.86236", "5.01584", "10.82974")
    for q, lower in zip(charges, floors, strict=True):
        stripped, _ = magnetic_matrix(q, {}, omit_clock_edges=True)
        inertia_lower_bound(stripped, edge_floor, arb(lower))
    print("无钟三角相位和删钟后的五荷恒秩一作用下界", triangle.arg().str(12), *floors)

    number = 8359
    height = acb.zeta_zero(number).imag
    clocks = geometry["phases"](height)
    model, *_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    layers = model["shell"]["shell_spectra"]((0, *charges), zero_index=number)
    center = sum(layers[0])/18
    lower = (".5283", "2.0860", "4.5925", "7.9183", "9.99999999")
    upper = (".57765", "2.28084", "5.02339", "8.67118", "10")
    counts = []
    for q, lo, hi in zip(charges, lower, upper, strict=True):
        lo, hi = arb(lo), arb(hi)
        if q != 6:
            full, trial = magnetic_matrix(q, clocks)
            inertia_lower_bound(full, edge_floor, lo)
            assert trial < hi < 9+edge_floor and lo < 9+edge_floor
        levels = chiral["ordered"]([v*abs(v-center) for v in layers[q]])
        below = sum(h < lo for h in levels)
        assert below + sum(h > hi for h in levels) == 19
        counts.append(below)
        print("绝对荷、实投影极小值包围和两步质量轻模数", q, lo.str(10), hi.str(10), below)
    assert counts == [3, 2, 3, 3, 3]


if __name__ == "__main__":
    verify()
