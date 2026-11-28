"""固定一个发射电流时，第120、163、177枚零点共用的最小接收站集。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


SITES = (0, 2, 3, 5, 6, 9, 11, 13)
OTHER = tuple(i for i in SITES if i != 6)


def area(z, w):
    return (z.conjugate()*w).imag


def check_witnesses(greens):
    for (n, g), extra in zip(greens.items(), (0, 2, 9), strict=True):
        for third in (0, extra):
            x, y, z = (g[6,j] for j in (third, 3, 5))
            assert max(area(y,z), area(z,x), area(x,y)) < arb("-1e-6")
        for j in OTHER:
            for i in SITES:
                assert g[i,j].real + {3:-1, 5:1}.get(j,0)*g[i,j].imag > arb(".00018")
        print("编号、可与站 3、5 搭配的第三站", n, extra)


def check_branches(greens):
    for g in greens.values():
        for i in SITES:
            for j in SITES:
                special = {i,j} in ({6,3}, {6,5})
                assert g[i,j].real > (arb("-.00085") if special else arb(".006"))
    assert arb(".006") / arb(".00685") > arb(".87")

    # 小 a_6：两枚编号各自除第六列外，投影都严格为正。
    def project(z, sign):
        return arb(".32")*z.real + sign*z.imag

    for index, sign in ((120,-1), (177,1)):
        g = greens[index]
        for i in OTHER:
            for j in OTHER:
                for z in (g[i,j], arb(".13")*g[6,j]+arb(".87")*g[i,j]):
                    assert project(z, sign) > arb(".00033")
    for i in SITES:
        assert project(greens[120][i,6],-1)+project(greens[177][i,6],1) > arb(".000398")

    # 大 a_6：每一枚的接收端都须在站 3 与站 5 分别有严格正权。
    for g in greens.values():
        for omitted, slope in ((3,5), (5,-5)):
            for i in SITES:
                for j in SITES:
                    if j != omitted:
                        z = arb(".87")*g[6,j] + arb(".13")*g[i,j]
                        assert z.real + slope*z.imag > arb(".002")


def coefficients(green):
    constant = area(green[6,3], green[6,5])
    assert constant < 0
    scale = -constant
    linear = {i:(area(green[i,3],green[6,5])+area(green[6,3],green[i,5]))/scale
              for i in OTHER}
    quadratic = {(i,j):(area(green[i,3],green[j,5])+area(green[j,3],green[i,5]))/(2*scale)
                 for i in OTHER for j in OTHER}
    return linear, quadratic


def exclude_two(greens):
    c, first = coefficients(greens[120])
    e, second = coefficients(greens[163])
    diff = {i:e[i]-c[i] for i in OTHER}
    assert min(c.values()) > 11 and min(diff.values()) > 1
    for i in OTHER:
        for j in OTHER:
            loss = arb("2.7")*(c[i]+c[j])/2
            difference_loss = arb(".24")*(c[i]*diff[j]+c[j]*diff[i])/2
            assert first[i,j]+loss > arb("2.8")
            assert second[i,j]-first[i,j]+difference_loss > arb(".14")
    assert arb(".13")/arb(".87") < arb(".15")
    assert 1/(1-arb("2.7")*arb(".15")) < arb("1.7")
    assert 1-arb(".24")*arb("1.7") > 0


def verify(replay):
    model = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 192
    greens = {n:model["laplacian"](5, acb.zeta_zero(n).imag).inv() for n in (120,163,177)}
    check_witnesses(greens)
    check_branches(greens)
    exclude_two(greens)
    print("一个固定发射端的接收站并集至少三站，各次支撑的交集至少两站，两界都能达到。")


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
