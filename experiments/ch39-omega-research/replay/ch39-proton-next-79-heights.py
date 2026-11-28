"""逐一检验第 178 至 256 枚正零点的八站非负静态相消。"""
from cmath import phase, rect
from itertools import combinations, combinations_with_replacement
from math import pi
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def separated(points):
    # 浮点只挑选分离方向。整张开半平面证书由随后严格区间内积验收。
    angles = sorted(phase(complex(z)) for z in points)
    gaps = [y-x for x,y in zip(angles, angles[1:]+[angles[0]+2*pi], strict=True)]
    index = max(range(len(gaps)), key=gaps.__getitem__)
    direction = rect(1, angles[index] + gaps[index]/2 + pi)
    normal = acb(round(10**7*direction.real), round(10**7*direction.imag))
    return all((normal.conjugate()*z).real > 0 for z in points)


def exclude_square(p, q, r, s):
    cells = [(0, 0)]
    while cells:
        numerator, depth = cells.pop()
        left, right = [arb(i)/2**depth for i in (numerator, numerator+1)]
        corners = [(1-t)*a+t*b for t in (left, right) for a,b in ((p,r), (q,s))]
        if not separated(corners):
            assert depth < 16, (p,q,r,s)
            cells.extend([(2*numerator, depth+1), (2*numerator+1, depth+1)])


def verify():
    root = Path(__file__).resolve().parent
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 192
    sites = (0, 2, 3, 5, 6, 9, 11, 13)
    easy, squares, feasible = [], [], []
    for index in range(178, 257):
        green = shell["laplacian"](5, acb.zeta_zero(index).imag).inv()
        if all(green[i,j].real > 0 for i in sites for j in sites):
            easy.append(index)
            continue
        x, y, z = [green[6,j] for j in (0,3,5)]
        areas = [(u.conjugate()*v).imag for u,v in ((y,z), (z,x), (x,y))]
        if all(t > 0 for t in areas) or all(t < 0 for t in areas):
            feasible.append(index)
            continue
        for i in sites:
            assert separated([green[i,j] for j in sites])
        # 矩阵厄米，交换两端取复共轭，单行的检验也排除了全部单列。
        for (i, k), (j, m) in combinations_with_replacement(combinations(sites, 2), 2):
            exclude_square(green[i,j], green[i,m], green[k,j], green[k,m])
        squares.append(index)
    assert len(easy) == 75 and squares == [223, 236, 254] and feasible == [205]
    print("第 178 至 256 枚的实部排除数、有限支撑排除编号、可行编号",
          len(easy), squares, feasible)


if __name__ == "__main__":
    verify()
