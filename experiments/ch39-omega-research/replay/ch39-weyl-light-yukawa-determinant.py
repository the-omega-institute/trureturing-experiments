"""在荷负三 Weyl 球内核验其他手征荷的全部轻带可延拓，计算 Yukawa 行列式障碍。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    native = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    cone = run_path(str(root / "ch39-first-chern-edge.py"))
    cone["verify"](root)
    ctx.prec = 256
    center = arb("80.928425691033686017")
    weights = [arb(".01615054064436668415"), arb("11.16447955978410845")]
    # 含近似中心误差，1e-12 盒覆盖前述以精确交点为心的 1e-13 闭球。
    width = arb("1e-12")
    rates = [abs(native["geo"]["theta"](d)-native["geo"]["theta"](244))
             for _, _, d, _ in native["geo"]["edges"] if d in (724,804)]
    assert len(rates) == 12 and center > 80
    for charge in (1, -4, 2, 6):
        matrix = native["laplacian"](charge, center)
        for (i, j), weight in zip(cone["EDGES"], weights, strict=True):
            assert abs(matrix[i,j]).contains(1)
            delta = weight-1
            matrix[i,i] += delta
            matrix[j,j] += delta
            matrix[i,j] *= weight
            matrix[j,i] = matrix[i,j].conjugate()
        eig = matrix.eig(multiple=True)
        assert all(z.imag.contains(0) for z in eig)
        eig = sorted((z.real for z in eig), key=float)
        # 钟边相位导数、软边导数及两条单位磁边的权导数。
        shift = width*(abs(charge)*sum(rates)+1+4)
        gap = eig[3]-eig[2]-2*shift
        assert gap > arb(".01")
        if charge == 6:
            fourth_gap = eig[4]-eig[3]-2*shift
            assert fourth_gap > arb(".01")
            print("荷六第四、第五带之间的谱隙下界", fourth_gap.lower().str(15))
        print("荷、整个闭球的第三第四带隙下界", charge, gap.lower().str(15))
    print("带电轻子、Majorana 行列式线次数模三", (-1) % 3, (-2) % 3)


if __name__ == "__main__":
    verify()
