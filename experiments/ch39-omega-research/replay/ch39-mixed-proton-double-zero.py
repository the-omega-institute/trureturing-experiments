"""同一固定真空、同一绕数权和同一正反等分面混合可在两个指定点精确相消。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


FIXED = {
    1: ".024378240692", 3: ".020066094900", 5: ".005519462214",
    7: ".004773759048", 16: ".003527630636", 33: ".007992424541",
    35: ".281359414361", 36: ".045737144634", 67: ".003333497658",
    69: ".019104903577",
}
SOLVED = (104, 38, 81, 105, 14)
CLASSES = ((1, 0), (0, 1), (1, -1))
NUMBER = 2_306_888_100_085


def verify(root):
    def load(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    fixed = load("height-independent-flavor-vacuum")
    model = load("reference-flavor-global-winding")
    lib = load("mincut-global-winding")
    shell = load("three-ray-chiral-boundary")
    angles = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    model["conjugation_symmetry"](root)
    ctx.prec = 224
    tri = lib["terms"](shell)
    sizes = [len(tri[k]) for k in CLASSES]
    assert sizes == [42, 42, 12]
    assert [sum(soft for _, soft in tri[k]) for k in CLASSES] == [3, 3, 0]

    x = arb("-.01046")
    d = angles["reference_angles"](flavor)[0]
    y = x-d
    height = acb.zeta_zero(NUMBER).imag
    phi = shell["geo"]["phases"](height)
    actual = (phi[724].arg(), phi[804].arg())
    assert abs(actual[1]-y) < arb("9e-6")
    assert abs(actual[0]-actual[1]-d) < arb("2e-7")
    delta = arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])

    # 五个未知绕数权由归一化和跨高度共用的四个梯度分量精确定义。
    weights = {q: arb(value) for q, value in FIXED.items()}
    columns = {q: [arb(1), *fixed["gradient_columns"](lib, tri, q, x, y)]
               for q in (*FIXED, *SOLVED)}
    mat = arb_mat([columns[q] for q in SOLVED]).transpose()
    assert not mat.det().contains(0)
    rhs = [int(j == 0)-sum(w*columns[q][j] for q, w in weights.items())
           for j in range(5)]
    weights.update(zip(SOLVED, mat.solve(arb_mat([[v] for v in rhs])).entries(), strict=True))
    assert min(weights.values()) > 0 and max(weights) == 105
    print("正绕数权", [(q, w.str(15)) for q, w in sorted(weights.items())])
    for j in range(5):
        residual = sum(w*columns[q][j] for q, w in weights.items())-int(j == 0)
        assert residual.contains(0) and residual.rad() < arb("1e-45")

    # 用该编号的真实软边亏额同时计算两点。正反等分把复因子投影到实部。
    coeffs = {q: lib["coeff"](tri, q, delta) for q in weights}

    def factors(point):
        a, b = point
        return [sum(w*(coeffs[q][j]*acb(0, q*t).exp()).real
                    for q, w in weights.items())/sizes[j]
                for j, t in enumerate((a, b, a-b))]

    reference = factors((x, y))
    numbered = factors(actual)
    differences = [(n-r)*10**6 for r, n in zip(reference, numbered)]
    mixmat = arb_mat([[1, 1, 1], reference, differences])
    det = mixmat.det()
    assert not det.contains(0)
    shares = mixmat.solve(arb_mat([[1], [0], [0]])).entries()
    for s, lo, hi in zip(shares, (".0766", ".0151", ".9082"), (".0767", ".0152", ".9083")):
        assert arb(lo) < s < arb(hi)
    for row, expected in (([1, 1, 1], 1), (reference, 0), (numbered, 0)):
        residual = sum(s*a for s, a in zip(shares, row))-expected
        assert residual.contains(0) and residual.rad() < arb("1e-36")
    slopes = [sum(w*(acb(0,q)*coeffs[q][j]*acb(0,q*t).exp()).real
                  for q,w in weights.items())/sizes[j] for j,t in enumerate((x,y,d))]
    gx = shares[0]*slopes[0] + shares[2]*slopes[2]
    gy = shares[1]*slopes[1] - shares[2]*slopes[2]
    direction = [actual[0]-x, actual[1]-y]
    length = sum(t*t for t in direction).sqrt()
    along = (gx*direction[0]+gy*direction[1])/length
    across = (gy*direction[0]-gx*direction[1])/length
    print("两钟梯度、双零连线与垂直方向导数",
          *(v.str(16) for v in (gx,gy,along,across)))
    assert abs(across) > arb(".1")
    # 两钟共同平移与改变钟差的窄窗，比较一个明确窗内点与精确零。
    shifted = factors((x+arb("1.9e-7"),y))
    witness = sum(a*b for a,b in zip(shares,shifted,strict=True))
    print("仅将第一钟转动1.9e-7弧度后的混合因子",witness.str(16))
    assert abs(witness) > arb("1e-8")
    print("参考真空三类实部", *(z.str(17) for z in reference))
    print("编号零点三类实部", *(z.str(17) for z in numbered))
    print("两处差乘以一百万", *(z.str(17) for z in differences))
    print("混合矩阵行列式与三类份额", *(z.str(17) for z in [det, *shares]))

    energy = sum(w*(lib["energy"](q, coeffs[q], *actual)
                    - lib["energy"](q, coeffs[q], x, y)) for q, w in weights.items())
    print("编号零点激发能", energy.str(18))
    assert 0 < energy < arb("1e-7")
    cap = arb(23)/25*arb(2)**-20/14**2
    assert 0 < delta < cap
    all_deficits = arb(0).union(cap)
    modes = [(q, [w*c for c in lib["coeff"](tri, q, all_deficits)]) for q, w in weights.items()]
    drift = model["global_cover"](lib, tri, weights, modes, x, y)
    assert drift < -3


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
