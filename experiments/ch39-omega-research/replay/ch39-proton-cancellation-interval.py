"""固定发射站时，认证三站及八站接收的荷五相消区间。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    endpoint = run_path(str(root / "ch39-proton-endpoint-speed.py"))
    ctx.prec = 192
    height = acb.zeta_zero(177).imag
    left = ("-.0417115028460", "-.0417115028458")
    right = (".1054886795107", ".1054886795110")
    sites = (3, 5, 13)
    area = endpoint["area"]

    def triangle(green):
        x, y, z = [green[6, j] for j in sites]
        return area(y, z), area(z, x), area(x, y)

    # 第一、第二个面积分别严格下降、上升。第三个面积在全区间严格为负。
    # 用相邻有理端点的区间包络覆盖整段高度，以包含矩阵求逆的误差。
    for k in range(-42, 106):
        for j in range(8):
            numerator = 8*k + j
            offset = (arb(numerator)/8000).union(arb(numerator + 1)/8000)
            green, slope = endpoint["green_and_slope"](shell, height + offset)
            x, y, z = [green[6, i] for i in sites]
            dx, dy, dz = [slope[6, i] for i in sites]
            assert area(dy, z) + area(y, dz) < 0
            assert area(dz, x) + area(z, dx) > 0
            assert area(x, y) < 0

    for index, endpoints, expected_signs in [(0, left, (1, -1)), (1, right, (-1, 1))]:
        for offset, sign in zip(endpoints, expected_signs, strict=True):
            green = shell["laplacian"](5, height + arb(offset)).inv()
            assert sign * triangle(green)[index] > 0
    print("包含第 177 枚正零点的相消分支，两端相对零点的高度位移分别在", left, right)

    # 在原三角形的两端各接上一段。r 是支持直线的方向，k 为穿过该直线的站。
    sites8 = (0, 2, 3, 5, 6, 9, 11, 13)
    sides = (
        (-47, -41, 0, 5, 3, 1, ("-.0459972646597", "-.0459972646594")),
        (105, 107, 2, 3, 5, -1, (".1064730230790", ".1064730230793")),
    )
    for start, stop, r, k, companion, inside_sign, bracket in sides:
        for numerator in range(32*start, 32*stop):
            offset = (arb(numerator)/32000).union(arb(numerator + 1)/32000)
            green, slope = endpoint["green_and_slope"](shell, height + offset)
            row = {i: green[6, i] for i in sites8}
            assert abs(row[r]) > 0
            for i in sites8:
                if i not in (r, k):
                    assert inside_sign * area(row[i], row[r]) > 0
            assert inside_sign * area(row[k], row[companion]) > 0
            assert area(slope[6, r], row[k]) + area(row[r], slope[6, k]) > 0
        for offset, sign in zip(bracket, (-1, 1), strict=True):
            green = shell["laplacian"](5, height + arb(offset)).inv()
            assert sign * area(green[6, r], green[6, k]) > 0
        print("八站接收时的边界站对与高度位移", (r, k), bracket)


if __name__ == "__main__":
    verify()
