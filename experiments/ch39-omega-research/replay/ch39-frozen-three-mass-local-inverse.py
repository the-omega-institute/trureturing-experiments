"""冻结 K、标尺后，实际中央附近的三个归一化质量共同反演两角及径向亏额。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


CENTERS = {-4: 2345666, 2: 964890}


def mass_map(model, charge):
    """返回实际中央亏额与 evaluate(724 弦角位移, 804 弦角位移, 绝对亏额)。

    evaluate 依次返回三个按参考半宽归一的质量及它们的 3×3 雅可比。
    两个角坐标均已乘本部门的荷，均以该部门实际零点的相位为原点。
    """
    ctx.prec = 384
    shell, pump, adj = model["shell"], model["pump"], model["adj"]
    height = acb.zeta_zero(CENTERS[charge]).imag
    phases = shell["geo"]["phases"](height)
    offset = shell["soft"]["point"]["study"]["offset"]
    actual = (arb(23)/25) * arb(2)**-20 / (height**2 + offset)
    base, edges = pump["clock_matrices"](shell, charge)
    reference = [edge*phases[d]**charge for edge, d in zip(edges, (724,804))]
    scale = arb("7.06539940386823251") if charge == -4 else arb(1)
    widths = model["tol"][:3] if charge == -4 else model["tol"][3:6]

    def evaluate(alpha, beta, deficit):
        terms = [edge*acb(0, angle).exp() for edge, angle in zip(reference, (alpha, beta))]
        L = base + sum((z + adj(z) for z in terms), 0*base) - deficit*model["radial"]
        directions = [acb(0,1)*(z-adj(z)) for z in terms] + [-model["radial"]]
        gram = acb_mat(3,3)
        derivatives = [acb_mat(3,3) for _ in directions]
        # 每个字符的最低单根只贡献 f(e)^2 P；约化预解式给出 dP。
        for frame in model["frames"]:
            levels, basis = model["eig"](adj(frame)*L*frame)
            cols = [acb_mat([[basis[i,j]] for i in range(basis.nrows())])
                    for j in range(basis.ncols())]
            projectors = [v*adj(v) for v in cols]
            e = levels[0]
            f = 1/(1 + (e/9)**512)
            df = -(arb(512)/9)*(e/9)**511*f**2
            S = adj(frame)*model["sources"]*model["K"]*scale
            gram += adj(S)*(f**2*projectors[0])*S
            reduced = sum((projectors[j]*(1/(e-levels[j])) for j in range(1,len(levels))),
                          0*projectors[0])
            for k, op in enumerate(directions):
                X = adj(frame)*op*frame
                de = (adj(cols[0])*X*cols[0])[0,0].real
                half = reduced*X*projectors[0]
                derivatives[k] += adj(S)*(2*f*df*de*projectors[0] + f**2*(half+adj(half)))*S
        levels, basis = model["eig"](gram)
        masses = [lev.sqrt() for lev in levels]
        J = arb_mat([[(adj(basis)*op*basis)[i,i].real/(2*masses[i]*widths[i])
                      for op in derivatives] for i in range(3)])
        return [mass/width for mass,width in zip(masses,widths)], J

    return actual, evaluate


def check(replay):
    model = run_path(str(replay/"ch39-two-zero-quark-window.py"))
    for q in CENTERS:
        actual, evaluate = mass_map(model, q)
        zero = arb(0)
        radius = arb("1e-12")
        angles = arb(0,radius.upper())
        loss = zero.union(arb("2e-14"))
        assert 0 < actual < arb("2e-14")
        _, center = evaluate(zero, zero, actual)
        _, box = evaluate(angles, angles, loss)
        identity = arb_mat([[int(i==j) for j in range(3)] for i in range(3)])
        error = center.inv()*box - identity
        norm = max(sum(abs(error[i,j]).upper() for j in range(3)) for i in range(3))
        assert norm < arb("4.40e-6" if q == -4 else "6.06e-5")
        # 任意两点间的平均雅可比也逐项落在 box 中。乘中央的逆后，
        # 它与恒等映射相距小于 1，因此三输出映射在整个凸盒单射。
        inverse = box.inv()
        radial_gain = 1 / sum(abs(inverse[2,j]) for j in range(3)).upper()
        discrepancy = actual * radial_gain
        assert radial_gain > arb("38.96133" if q == -4 else "35.73931")
        assert discrepancy > arb("1.99100e-17" if q == -4 else "9.28911e-17")
        minors = [box[i,0]*box[j,1]-box[j,0]*box[i,1] for i,j in combinations(range(3),2)]
        assert all(minor > 0 for minor in minors)
        print(q, "实际亏额",actual.str(15),"中央行列式",center.det().str(12),
              "预条件误差上界",norm.str(14),"径向增益下界",radial_gain.str(14),
              "平坦质量差下界",discrepancy.str(14))
        print(q, "角子式 (01,02,12)",[minor.str(14) for minor in minors])


if __name__ == "__main__":
    check(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
