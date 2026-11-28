"""任意固定的八站非负端点电流在前三枚可相消零点中至多消去一枚。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 192
    sites = (0, 2, 3, 5, 6, 9, 11, 13)
    greens = {}
    for index in (120, 163, 177):
        green = shell["laplacian"](5, acb.zeta_zero(index).imag).inv()
        for i in sites:
            for j in sites:
                exceptional = {i, j} in ({6, 3}, {6, 5})
                assert green[i, j].real > (arb("-.00085") if exceptional else arb(".006"))
        greens[index] = green

    # 若三种配对中任意一对都很小，取共轭使 p>=q 后即有 a_6>0.8。
    c = (arb(".006") - arb(".00024")) / arb(".00685")
    assert c - (1-c)**2/(2*c) > arb(".8")
    cases = (
        (120, 163, (-9304,5148,10000,-4990), ".000015"),
        (120, 177, (845,9060,-355,-10000), ".00024"),
        (163, 177, (1717,8931,-831,-10000), ".00024"),
    )
    for first, second, weights, minimum in cases:
        x, y, u, v = [arb(k)/10000 for k in weights]
        norm_sum = (x*x + y*y).sqrt() + (u*u + v*v).sqrt()
        for i in sites:
            for j in sites:
                a, b = [(4*g[6, j] + g[i, j])/5 for g in (greens[first], greens[second])]
                projection = x*a.real + y*a.imag + u*b.real + v*b.imag
                assert projection > arb(minimum)*norm_sum
        print("编号对及两个振幅模长最大值的严格下界", first, second, minimum)


if __name__ == "__main__":
    verify()
