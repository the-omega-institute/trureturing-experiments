"""三窗零点上的全温度带荷壳热迹相对首枚的上界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


root = Path(__file__).parent
thermal = run_path(str(root / "ch39-shell-thermal-time-weights.py"))
shell = thermal["shell"]


def mass_floors(charge, height):
    edges = shell["geo"]["edges"]
    clock_vertices = {6, 7, 8}
    assert all(i in clock_vertices or j in clock_vertices
               for i, j, d, _ in edges if d in (724, 804))
    core = [i for i in range(19) if i not in clock_vertices]
    full = shell["laplacian"](charge, height)
    principal = acb_mat([[full[i,j] for j in core] for i in core])
    eigen = principal.eig(multiple=True)
    assert len(eigen) == 16 and all(v.imag.contains(0) for v in eigen)
    lower = sorted(v.real.lower() for v in eigen)
    upper = sorted(v.real.upper() for v in eigen)

    offset = shell["soft"]["point"]["study"]["offset"]
    assert offset >= 0
    change = arb(23)/25*arb(2)**-20/(height**2 + offset)
    assert 6*change < arb("1e-9")
    # 135 条无向边中三条辐边减重，故中心为 15−δ/3。
    pole, tuned = shell["soft"]["pole"], shell["soft"]["star"]["orbit"]
    assert len(edges) == 135
    assert sum(pole in (i,j) and (i in tuned or j in tuned) for i,j,_,_ in edges) == 3
    values = []
    for j in range(19):
        lo = max(arb(0), (lower[j-3]-6*change).lower()) if j >= 3 else arb(0)
        hi = min(arb(36), (upper[j]+6*change).upper()) if j < 16 else arb(36)
        assert lo <= hi
        if hi < 15 or lo > 15:
            floor = min(abs(lo*(lo-15)).lower(), abs(hi*(hi-15)).lower())
            values.append(max(arb(0), (floor-12*change).lower()))
        else:
            values.append(arb(0))
    return [max(arb(9), v) for v in sorted(values)[3:]]


def verify():
    ctx.prec = 224
    height = acb.zeta_zero(8359).imag
    reference = thermal["heavy_spectrum"](8359)
    for charge, masses in reference.items():
        lower = mass_floors(charge, height)
        assert len(lower) == len(masses) == 16
        for power, bound in ((1, "1.680"), (2, "2.147")):
            thermal["below_ratio"]([m**power for m in masses], [m**power for m in lower], arb(bound))
    assert arb("1.680")**36 < 195557000 and arb("2.147")**24 < 195557000
    print("全部三窗候选相对首枚的单热迹倍率上界", 1.680, 2.147)
    print("质量一次、允许质量平方时已排除的因子数上限", 36, 24)


if __name__ == "__main__":
    verify()
