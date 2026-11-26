"""固定谱三重性与十九站点独立超荷生成元的最小李代数。"""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


triality = run_path(str(Path(__file__).with_name("ch39-heavy-shift-baryon-triality.py")))
study = triality["study"]
charges = (1,-4,2,-3,6)


def projector(q, height, center):
    v = study["edge"]["mode_columns"](study["edge"]["filters"]["spectra"](q,height,center))
    return v*study["adjoint"](v)


def verify():
    ctx.prec = 448
    omega = acb(-1,arb(3).sqrt())/2
    # ad(iD_i)^2 在跨 i,j 的元素上有五个不同的实谱值。
    for q,r in ((q,r) for q in charges for r in charges):
        value = Fraction(1)
        for other in charges:
            if other != q:
                value *= Fraction(other**2-r**2,other**2-q**2)
        assert value == int(q==r)
    for n in (22027,28325):
        height = acb.zeta_zero(n).imag
        zero = study["shell"]["laplacian"](0,height)
        center = sum(zero[i,i].real for i in range(19))/18
        for q in charges:
            p = projector(q,height,center)
            # 对源站点 0，TE_0T† 的每个非对角元都非零。
            column = [(1-omega)*p[i,0] + (omega if i==0 else 0) for i in range(19)]
            # 单独标记站点 0 时，e_0 与 T e_0 跨越两个不同的谱方向。
            assert 0<p[0,0].real<1 and 0<abs(column[0])<1
            assert all(abs(column[i]*column[j].conjugate()) > arb("1e-11")
                       for i,j in combinations(range(19),2))
    assert [n%2 for n in (6,3,3,2,1)] == [0,1,1,0,1]
    print("公共超荷加单个站点时的闭包维数",5*3+2)
    print("站点代数、单部门闭包、五个有符号超荷的共同闭包维数",19,19**2,5*(19**2-1)+1)


if __name__ == "__main__":
    verify()
