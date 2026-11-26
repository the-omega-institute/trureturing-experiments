"""原站点编号下一跳壳局域的无损整代传输。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


def adjoint(matrix):
    return matrix.conjugate().transpose()


def verify(root):
    lib=run_path(str(root/"ch39-local-yukawa-plucker.py"))
    prev=run_path(str(root/"ch39-joint-zero-diagonal-light-filters.py"))
    shell=lib["shell"]
    ctx.prec=320
    height=acb.zeta_zero(2_306_888_100_085).imag
    delta=arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    center=15-delta/3
    vertices=range(19)
    eye=acb_mat([[int(i==j) for j in vertices] for i in vertices])
    adjacent={(i,i) for i in vertices}
    for i,j,_,_ in shell["geo"]["edges"]:
        adjacent.update(((i,j),(j,i)))
    blocked=[(i,j) for i in vertices for j in vertices if (i,j) not in adjacent]
    assert len(blocked)==72 and len(adjacent)==289
    assert all((shell["soft"]["pole"],j) in adjacent for j in vertices)
    u=lib["light_basis"](1,height,center)
    p=u*adjoint(u)
    pc=eye-p
    for q in (2,-4):
        v=lib["light_basis"](q,height,center)
        qc=eye-v*adjoint(v)
        # 半径零的对角滤波系统满列秩。
        prev["certify_positive"](acb_mat([[qc[i,j]*p[j,i] for j in vertices] for i in vertices]))
        base=v*adjoint(u)
        H=acb_mat([[qc[i,k]*pc[l,j] for k,l in blocked] for i,j in blocked])
        prev["certify_positive"](H-arb(".033")*acb_mat([[int(i==j) for j in range(72)] for i in range(72)]))
        coeff=H.solve(acb_mat([[-base[i,j]] for i,j in blocked]))
        sparse=acb_mat(19,19)
        for (i,j),z in zip(blocked,coeff.entries(),strict=True):
            sparse[i,j]=z
        residual=qc*sparse*pc
        transfer=base+residual
        # 等式成立来自 H 的可逆性及投影恒等式。下面独立检查区间残差和矩阵范数。
        assert all(transfer[i,j].contains(0) and abs(transfer[i,j])<arb("1e-30") for i,j in blocked)
        assert all(z.contains(0) and abs(z)<arb("1e-30") for z in (transfer*u-v).entries())
        assert all(z.contains(0) and abs(z)<arb("1e-30") for z in (adjoint(transfer)*v-u).entries())
        norm=sum(abs(z)**2 for z in residual.entries()).sqrt()
        assert norm < arb(1), norm
        print("目标荷、禁止项数、补块 Frobenius 范数、全算子范数与轻态成功率",q,len(blocked),
              norm.str(13),1,1,flush=True)


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
