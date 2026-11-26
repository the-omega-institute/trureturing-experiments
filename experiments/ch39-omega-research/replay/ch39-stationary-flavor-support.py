"""局部密度势的空站点场方程与顶底、粲顶层级的必要支撑。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def positive_definite(matrix):
    # Hermitian Sylvester criterion; a strictly signed pivot is required.
    for size in range(1,matrix.nrows()+1):
        leading = acb_mat([[matrix[i,j] for j in range(size)] for i in range(size)])
        det = leading.det()
        assert det.imag.contains(0) and not det.real.contains(0)
        if det.real < 0:
            return False
    return True


def verify():
    root = Path(__file__).parent
    flavor = run_path(str(root/"ch39-minimum-flavor-support.py"))
    cones = run_path(str(root/"ch39-ten-site-charm-exclusion.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 192
    up,down = flavor["flavor_maps"](shell,11431)
    lap = shell["laplacian"](3,acb.zeta_zero(11431).imag)
    sites = shell["geo"]["sites"]
    rot = [sites.index((y,z,x)) for x,y,z in sites]
    ref = arb(".8573")/arb(".01410")
    adjoint = flavor["adjoint"]
    flavor_form = adjoint(down)*down-3*adjoint(up)*up/ref**2
    for size in (13,14,15):
        weight = arb(".01") if size == 13 else arb(1)
        margin = arb(".00000036") if size == 13 else arb("1e-10")
        target = arb(".40") if size == 14 else arb(".015")
        groups = cones_used = 0
        for group in combinations(range(19),size):
            if group != min(flavor["orbit"](group,rot)):
                continue
            groups += 1
            border = acb_mat([[lap[i,j] for j in group] for i in range(19) if i not in group])
            constrained = flavor["submatrix"](flavor_form,group)+weight*adjoint(border)*border
            for j in range(size):
                constrained[j,j] -= margin
            if size != 15 and positive_definite(constrained):
                continue
            assert size != 13
            us = acb_mat([[up[i,j] for j in group] for i in range(9)])
            cones["check_group"](constrained,us,target)
            cones_used += 1
        assert groups == {13:9054,14:3876,15:1296}[size]
        print("站点数、循环代表数、需用粲质量约束的代表数、空站残差下界",
              size,groups,cones_used,(margin/weight).sqrt().str(8),flush=True)


if __name__ == "__main__":
    verify()
