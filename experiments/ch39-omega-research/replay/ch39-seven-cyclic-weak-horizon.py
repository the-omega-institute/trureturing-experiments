"""固定第五壳源后，任何同一个三循环等变酉输运都不能在第七至十三壳逐壳消去弱差。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


CANDIDATES = ((0,9,15), (0,9,15), (0,15), (0,15), (0,9,15), (0,9,15), (0,9,15))


def adj(a):
    return a.conjugate().transpose()


def norm(a):
    return sum(abs(z)**2 for z in a.entries()).sqrt()


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    ninth = lib("ninth-weak-two-line-threshold")
    squared, _, _, _, cover = ninth["verify"](root)
    native = lib("native-rh-ray-pvm")
    graph = lib("three-ray-chiral-boundary")
    characters = lib("cyclic-two-characters-nine")["characters"]
    trivial = characters(graph["geo"]["sites"])[0][0]
    assert trivial.ncols() == 7
    gram = native["gram_cholesky"](root)
    floor, heavy, top = levels
    ctx.prec = 256
    assert squared[0].sqrt().overlaps(floor)
    planck = arb("1.22089e19")
    mz = arb("91.1876")
    q = arb(3000)
    far = arb(".000020")
    near = arb(".020483945")
    lower = arb("-1.4e-11")
    reference = arb("-.074256577194663365545")
    cover((q/mz)**2, far)
    cover(floor**2/(mz/planck)**2, near)
    weighted = (fifth*columns)*acb_mat([[(heavy-floor).sqrt(), 0], [0, (top-floor).sqrt()]])
    required = norm(adj(trivial)*weighted)
    allowed = (q/planck-floor).sqrt()
    all_rows = []
    upper_weak = []
    for j, ids in enumerate(CANDIDATES):
        _, basis = lib("adjacent-pvm-transport")["spectrum"](graph, 6, 2_306_888_100_092+j)
        weights = [z.real for z in native["readout"](basis, range(19), [14,8,7,6], gram[3], 3)]
        assert all(not v.contains(0) for v in weights) and sum(weights).contains(1)
        assert all(weights[i] > 0 for i in ids)
        upper = reference + sum(w*(far if i in ids else near if w>0 else lower)
                                for i,w in enumerate(weights))
        assert upper<0
        upper_weak.append(upper.upper())
        left = adj(basis)
        rows = {}
        for i in ids:
            row = acb_mat([[left[i,k] for k in range(19)]])
            b = row*trivial
            e = sum(abs(z).upper()**2 for z in (row-b*adj(trivial)).entries()).sqrt()
            rows[i] = (b, (allowed+e*norm(weighted)).upper())
        all_rows.append(rows)
    # 三循环等变输运保持 Z 的 Frobenius 范数。用逆矩阵的列和 BZ 的行
    # 将 Z 写成七个秩一外积之和，再应用三角不等式。
    upper_norm = []
    worst = None
    for ids in product(*CANDIDATES):
        data = [rows[i] for rows,i in zip(all_rows, ids, strict=True)]
        matrix = acb_mat([[b[0,k] for k in range(7)] for b,_ in data])
        inv = matrix.inv()
        bound = sum(sum(abs(inv[k,j])**2 for k in range(7)).sqrt()*cap
                    for j,(_,cap) in enumerate(data)).upper()
        if worst is None or bound > worst[0]:
            worst = (bound, ids)
        upper_norm.append(bound)
    assert len(upper_norm)==972
    print("七壳候选全高时的弱差统一上界、枚举数、平凡字符源范数、低线所能容许的最大范数",
          max(upper_weak).str(14),len(upper_norm),required.lower().str(15),
          max(upper_norm).str(15), worst[1], flush=True)
    assert all(x.is_finite() for x in upper_norm)
    assert max(upper_norm)<arb(".626")<arb(".7")<required
    # 对不要求等变的 U，该差同时界住 ||[EE*,U]|| 与到所有等变酉矩阵的距离。
    separation = (required-max(upper_norm))/norm(weighted)
    print("七壳弱解到三循环等变酉族的谱范数距离下界", separation.lower().str(14), flush=True)
    assert separation>arb(".07")


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))
