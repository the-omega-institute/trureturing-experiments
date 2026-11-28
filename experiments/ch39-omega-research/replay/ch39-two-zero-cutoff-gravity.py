"""在两个相邻实际高度复用归一有符号截止，同时匹配牛顿与真空读数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx, fmpz_mat
import sympy as sp


def neutral_levels(shell, delta):
    """用有理软边的整数特征多项式及 Weyl 界包住十八条非零谱线。"""
    scale = 2**192
    loss = (scale*delta).mid().floor().unique_fmpz()
    matrix = fmpz_mat(19, 19)
    tuned = shell["soft"]["star"]["orbit"]
    pole = shell["soft"]["pole"]
    degree = 0
    for i, j, _, _ in shell["geo"]["edges"]:
        spoke = pole in (i, j) and (i in tuned or j in tuned)
        degree += int(spoke)
        weight = scale-loss if spoke else scale
        matrix[i, i] += weight
        matrix[j, j] += weight
        matrix[i, j] -= weight
        matrix[j, i] -= weight
    assert degree == 3
    # 三边星的拉普拉斯算子范数是 4。常数向量在原矩阵和近似矩阵中均给出零模。
    error = 4*abs(delta-arb(loss)/scale)
    roots = matrix.charpoly().complex_roots()
    vals = []
    zero_count = 0
    for value, multiplicity in roots:
        assert value.imag.contains(0)
        if value.real.is_zero():
            zero_count += multiplicity
        else:
            level = value.real/scale + arb(0, error.upper())
            assert level > 0
            vals.extend([level]*multiplicity)
    assert zero_count == 1 and len(vals) == 18
    return vals


def gravity_fields(root, numbers=None, *, heights=None, tight_neutral=False):
    single = run_path(str(root / "ch39-signed-cutoff-gravity.py"))
    lib = run_path(str(root / "ch39-chiral-minimal-neutral-gravity.py"))
    chiral = lib["chiral"]
    model, *_ = chiral["source"]["mean"]["inputs"]()
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 320
    flat = sp.zeros(19)
    for i, j, _, _ in shell["geo"]["edges"]:
        flat[i, i] += 1
        flat[j, j] += 1
        flat[i, j] -= 1
        flat[j, i] -= 1
    intervals = flat.charpoly().as_poly().intervals(eps=sp.Rational(1, 10**50))
    neutral = [arb(str(lo)).union(arb(str(hi)))
               for (lo, hi), multiplicity in intervals for _ in range(multiplicity)]
    assert len(neutral) == 19 and neutral[0].is_zero()
    def levels(h):
        delta = arb(23)/25 * arb(2)**-20 / (h*h+shell["soft"]["point"]["study"]["offset"])
        c = 15-delta/3
        bands = {}
        for q in (0, *model["charges"]):
            if q == 0:
                vals = (neutral_levels(shell, delta) if tight_neutral else
                        [v-arb(0).union(4*delta) for v in neutral[1:]])
            else:
                eig = shell["laplacian"](q, h).eig(multiple=True)
                vals = [z.real for z in eig]
                assert len(vals) == 19 and all(z.imag.contains(0) for z in eig)
            bands[q] = chiral["ordered"]([v*abs(v-c) for v in vals])
        rows = lib["fields"](bands)
        assert len(rows) == 99 and all(z > 0 for _, _, z in rows[1:])
        return rows
    if numbers is None:
        n = single["NUMBER"]
        numbers = (n, n + 1)
    if heights is None:
        heights = [acb.zeta_zero(n).imag for n in numbers]
    return single, lib, [levels(h) for h in heights]


def adjacent_fields(root):
    single, lib, rows = gravity_fields(root)
    return single, lib, *rows


def verify():
    root = Path(__file__).resolve().parent
    single, lib, rows0, rows1 = adjacent_fields(root)
    target = arb("1e-122")
    def read(rows, x):
        return single["moments"](rows, lib["expint"], x)
    _, sites, weights = single["capped"](rows0, lib["expint"], target, return_measure=True)
    mixed = [sum(p * v for p, v in zip(weights, axis))
             for axis in zip(*(read(rows1, x) for x in sites))]
    assert arb("1.34998e-7") < mixed[0]-1 < arb("1.35001e-7")
    assert arb("2.27669e-8") < mixed[1]/mixed[0]**2 < arb("2.27672e-8")
    print("直接复用单高度最优解，下一枚的牛顿倒数相对偏差及无量纲曲率",
          (mixed[0]-1).str(16), (mixed[1]/mixed[0]**2).str(16))
    def column(x):
        k0, w0 = read(rows0, x)
        k1, w1 = read(rows1, x)
        return [1, k0, w0, k1-k0, w1-w0]
    # 五列可同时指定两个高度的 K,W；将两背景求差以保留目标上的共同输入。
    sites = [arb(x) for x in ("7.543120063e-10", "3.556480306e-6", ".02894266125", ".54", "1")]
    basis = arb_mat([column(x) for x in sites]).transpose()
    assert not basis.det().contains(0)
    weights = basis.solve(arb_mat([[1], [1], [target], [0], [0]])).entries()
    assert all(w > 0 or w < 0 for w in weights)
    cost = sum(-p for p in weights if p < 0)
    assert arb("338016.30") < cost < 338017
    print("五个共同截止的系数与负量", *[p.str(16) for p in weights], cost.str(16))


if __name__ == "__main__":
    verify()
