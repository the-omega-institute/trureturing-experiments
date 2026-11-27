"""一般同谱输运从第七壳起不能连续七壳满足真空项与零弱差。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def prepare(root, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    floor, _, _ = levels
    graph = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    frames, _ = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])
    candidates = lib("seven-cyclic-weak-horizon")["CANDIDATES"]
    gram = native["gram_cholesky"](root)
    numbers = range(2_306_888_100_092, 2_306_888_100_099)
    fields = lib("three-index-free-atoms")["broken_fields"](root, numbers)[2]
    _, _, _, _, cover = lib("ninth-weak-two-line-threshold")["verify"](root)
    ctx.prec = 256
    planck, mz = arb("1.22089e19"), arb("91.1876")
    near = arb(".020483945")
    far = arb(".000020")
    cover((floor*planck/mz)**2, near)
    cover((arb(3000)/mz)**2, far)
    limits = {700: arb(".000290"), 2000: arb(".000035")}
    for energy, upper in limits.items():
        cover((arb(energy)/mz)**2, upper)
    refined = (700, 700, 2000, 2000, 700, 700, 700)
    reference = arb("-.074256577194663365545")
    lower = arb("-1.4e-11")
    shells = []
    weak_upper = []
    for number, ids, energy, field in zip(numbers, candidates, refined, fields, strict=True):
        _, basis = lib("adjacent-pvm-transport")["spectrum"](graph, 6, number)
        rows = adj(basis)*frames[0]
        w = [z.real for z in native["readout"](basis, range(19), [14,8,7,6], gram[3], 3)]
        assert all(not v.contains(0) for v in w) and all(w[i] > 0 for i in ids)
        # 真实谱简单，磁拉普拉斯与坐标三循环对易；平凡字符中的谱行恰有七条。
        triv = [i for i, row in enumerate(rows.tolist()) if sum(abs(z)**2 for z in row) > arb(".9")]
        assert len(triv) == 7 and set(ids).issubset(triv)
        nontrivial = [i for i in range(19) if i not in triv]
        assert all(w[i] < 0 for i in nontrivial)
        # 一个候选低于 3000。若只有第一条能低于 3000，第一条还须低于 energy。
        for first in (far, limits[energy]):
            corrections = [first if i == 0 else far if i in ids else near if v > 0 else lower
                           for i, v in enumerate(w)]
            upper = reference + sum(v*t for v, t in zip(w, corrections, strict=True))
            assert upper < arb("-.001")
            weak_upper.append(upper.upper())
        pos = [i for i in triv if w[i] > 0]
        c = -sum(g for g, _, z in field)/(8*arb.pi())
        error = sum(abs(g)*z for g, _, z in field)/(4*arb.pi())
        assert c > 0 and all(z >= 0 for _, _, z in field)
        # 对 a>=0，error*a² <= (c/100)*a⁴ + 25*error²/c。
        allowance = sum(-v for v in w if v < 0)*25*error**2/c
        holder = sum((-w[i])**(-arb(1)/3) for i in nontrivial)**-3
        positive_rows = acb_mat([rows.tolist()[i] for i in pos])
        shells.append((rows, pos, positive_rows, w, c, error, allowance, holder))
    print("十四个弱核排除界中的最大值",max(weak_upper).str(15),flush=True)
    return candidates, shells, refined


def upper_bound(levels, shells, ids, caps):
    floor, first, top = levels
    planck = arb("1.22089e19")
    radii = [((arb(q)/planck-floor).sqrt()).upper() for q in caps]
    selected = acb_mat([s[0].tolist()[i] for s, i in zip(shells, ids, strict=True)])
    inv = selected.inv()
    trace_cap = sum(r*sum(abs(inv[k,j])**2 for k in range(7)).sqrt()
                    for j, r in enumerate(radii))**2
    # 对补空间十二行，截止之和不小于两列总迹减去平凡块迹上界。
    remaining = (top+first-2*floor-trace_cap).lower()
    assert remaining > 0
    bounds = []
    for _, pos, positive_rows, weights, c, error, allowance, holder in shells:
        positive = arb(0)
        for i, row in zip(pos, (positive_rows*inv).tolist(), strict=True):
            row_cap = sum(abs(z)*r for z, r in zip(row, radii, strict=True))**2
            amplitude = min(top.upper(), (floor+min(trace_cap.upper(), row_cap.upper())).upper())
            positive += weights[i]*(c*amplitude**4 + error*amplitude**2)
        negative_cost = arb(99)/100*c*holder*remaining**4
        bounds.append((positive+allowance-negative_cost).upper())
    assert bounds[3] < arb("-.00010")  # 第七壳起的第四个读数，即第十壳。
    return min(bounds)


def experiment(root, levels):
    candidates, shells, refined = prepare(root, levels)
    excluded = []
    first_lines = (0,)*7
    for ids in product(*candidates):
        if ids != first_lines:
            excluded.append(upper_bound(levels, shells, ids, (3000,)*7))
    assert len(excluded) == 971 and max(excluded) < arb("-.0024")
    final = upper_bound(levels, shells, first_lines, refined)
    assert final < arb("-.00010")
    print("971 种选线中的最大真空上界、全第一线细化后的上界",
          max(excluded).str(15), final.str(15), flush=True)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root, source["verify"](root, return_operator="spectral")[2])
