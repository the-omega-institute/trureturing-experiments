"""两荷独立的十九维酉运输将右邻训练谱权送至右外壳的精确夸克拟合。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx

CENTER = 2_306_888_100_085
GEOMETRY = (".007", ".045", ("-.19637446", "-.07341222", "-.84285237", ".49562850"))
# 各行给定十二个有理谱权的分子；其余七处由六个四站坐标和总和唯一确定。
SEEDS = (
    (None, 764867989137, 1922974962, 1922974962, 3211463956, 435428210279, None,
     None, 1922974962, 1831212075, 947547885715, 2097188069, 832204891889,
     947547885715, 1922974962, None, None, None, None),
    (None, 153766768075, 60206619486, 153766768075, 153766768075, None,
     21296524270, 153766768075, 153766768075, None, 153766768075, 17781895,
     293761665526, 15324547, 7345128, None, None, None, None),
)


def coordinates(h):
    p, d = h[0, 0].real, h[1, 1].real
    z, helicity = (h[2, 2].real+h[3, 3].real)/2, h[2, 3].imag
    return [p+d, 2*h[0, 1].real, -2*h[0, 1].imag, p-d, z-helicity, z+helicity]


def reconstruct(rows, goal, original, seed):
    rhs = [*coordinates(goal), sum(original)]
    cols = [[*coordinates(p), arb(1)] for p in rows]
    unknown = [i for i, n in enumerate(seed) if n is None]
    assert len(cols) == len(seed) == 19 and len(unknown) == 7
    y = [arb(n or 0)/10**12 for n in seed]
    rhs = [b-sum(col[j]*v for col, v in zip(cols, y, strict=True)) for j, b in enumerate(rhs)]
    mat = arb_mat([cols[i] for i in unknown]).transpose()
    assert not mat.det().contains(0)
    chosen = mat.solve(arb_mat([[b] for b in rhs]), algorithm="precond")
    for i, v in zip(unknown, chosen.entries(), strict=True):
        y[i] = v
    assert all(0 < v < 1 for v in y) and (sum(y)-sum(original)).contains(0)
    source = sorted((v.lower() for v in original), reverse=True)
    dest = sorted((v.upper() for v in y), reverse=True)
    margin = min((sum(source[:k])-sum(dest[:k])).lower() for k in range(1, 19))
    assert margin > 0
    raw = sum((p*v for p, v in zip(rows, y, strict=True)), acb_mat(4, 4))
    assert (raw-goal).contains(acb_mat(4, 4))
    return raw, margin


def verify(root):
    emb = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    raw = run_path(str(root / "ch39-free-four-site-native-positive-spectrum.py"))
    diag = run_path(str(root / "ch39-complex-four-site-unit-quark-amplitudes.py"))["diag"]
    ck = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))["standard_matrix"]
    projections = run_path(str(root / "ch39-three-gram-source-flavor.py"))["projections"]
    source = CENTER+1
    _, original, _ = emb["verify_nodes"](root, ((source, ".0294"),), ".1", (),
                                         ((6, 7, 10, 13), (6, 7, 10, 12)))
    sectors = raw["native"](root, CENTER+2)
    ctx.prec = 384
    adj = emb["adj"]
    mix = acb_mat(ck(*[arb(t).sin() for t in (".22704", ".04275", ".00377")], arb("1.139")))
    mass = [list(map(arb, s)) for s in ((".00000615", ".00311", ".8616"), (".0000135", ".000268", ".01401"))]
    a, b = [diag([y/m[-1] for y in m]) for m in mass]
    a = adj(mix)*a*mix
    basis = raw["target_basis"](a, b*arb(".1"))*diag([acb(0, 1), 1, acb(0, -1)])
    assert (adj(basis)*basis).contains(diag([1]*3))
    u, v, intrinsic, rotation = emb["common_embedding"](*GEOMETRY)
    frame = rotation*intrinsic
    assert (adj(frame)*frame).contains(diag([1]*3))
    proj = []
    for charge, rows, m, block, seed in zip((-4, 2), sectors, mass, (a, b), SEEDS, strict=True):
        target = adj(basis)*block*basis*m[-1]
        four = rotation*emb["extend_target"](target, u, v)*adj(rotation)
        physical, margin = reconstruct(rows, four, original[source, charge], seed)
        compressed = adj(frame)*physical*frame
        diff = compressed-target
        assert diff.contains(acb_mat(3, 3))
        error = max(abs(z).upper() for z in diff.entries())
        print("压缩矩阵区间误差", error.str(10), flush=True)
        assert error < arb("1e-38")
        levels, p = projections(compressed*compressed)
        assert all(z.overlaps(y*y) for z, y in zip(levels, m, strict=True))
        print("荷、十八个累计主序幅差的最小下界", charge, margin.str(15), flush=True)
        proj.append(p)
    up, down = proj
    assert all((p*d).trace().real.overlaps(abs(mix[i, j])**2)
               for i, p in enumerate(up) for j, d in enumerate(down))
    j = (up[0]*down[0]*up[1]*down[1]).trace().imag
    assert j > 0 and j.overlaps((mix[0, 0]*mix[1, 1]*mix[0, 1].conjugate()*mix[1, 0].conjugate()).imag)
    print("六质量、九个 CKM 模平方及带号 J 核验通过", j.str(15), flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
