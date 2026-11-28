"""枚举前 177 枚正零点，认证八站独立非负电流精确相消的全部编号。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 192
    sites = (0, 2, 3, 5, 6, 9, 11, 13)
    separating_normals = {36: (48, 88), 68: (71, -71), 83: (67, -75), 148: (59, -81)}
    simple = []
    witnesses = []
    for index in range(1, 178):
        height = acb.zeta_zero(index).imag
        green = shell["laplacian"](5, height).inv()
        if all(green[i, j].real > 0 for i in sites for j in sites):
            simple.append(index)
            continue
        if index in separating_normals:
            for i in sites:
                for j in sites:
                    exceptional = {i, j} in ({6, 3}, {6, 5})
                    lower = arb("-.00085") if exceptional else arb(".006")
                    assert green[i, j].real > lower
            # 实部相消时互换端点可取 a_6>0.86，此处只用 a_6>=0.8。
            x, y = separating_normals[index]
            for i in sites:
                for j in sites:
                    vertex = (4*green[6, j] + green[i, j])/5
                    assert (x*vertex.real + y*vertex.imag)/100 > arb(".00045")
            continue
        assert index in (120, 163, 177)
        x, y, z = [green[6, j] for j in (0, 3, 5)]
        areas = [(u.conjugate()*v).imag for u, v in ((y,z), (z,x), (x,y))]
        assert all(a > 0 for a in areas) or all(a < 0 for a in areas)
        witnesses.append(index)

    assert len(simple) == 170 and witnesses == [120, 163, 177]
    print("八站独立非负电流，实部排除数、半平面排除编号、精确相消编号",
          len(simple), sorted(separating_normals), witnesses)
    heights = [acb.zeta_zero(i).imag for i in (176, 177, 178)]
    assert heights[1]-heights[0] > arb("1.5456")
    assert heights[2]-heights[1] > arb("2.0419")


if __name__ == "__main__":
    verify()
