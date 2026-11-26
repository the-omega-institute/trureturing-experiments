"""联合近真空零点上，对角后选择与逐站保迹通道传递整个轻空间的界。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


ANGLES = {
    2: ("0", ".4321844991", ".2363907042", "2.2383811481", "-1.8650530527", "-2.9528038630", ".1886234598"),
    -4: ("0", "-.2887024361", ".1155735073", "1.0544141993", "-1.2682337844", "-.1157995278", "-.1103797503"),
}
EPSILON = arb(1) / 10**11
CONDITIONAL = {2: (".04609632",".04609634"), -4: (".05181978",".05181981")}
COMPLETE = {2: (".4674978685",".4674978688"), -4: (".4525132126",".4525132128")}


def certify_positive(matrix):
    """逐项作 Hermitian LDL* 消元。全部实枢轴有正下界。"""
    n=matrix.nrows()
    work=[[matrix[i,j] for j in range(n)] for i in range(n)]
    for k in range(n):
        pivot=work[k][k].real
        assert pivot>0
        for i in range(k+1,n):
            f=work[i][k]/pivot
            for j in range(k+1,i+1):
                work[i][j]-=f*work[j][k].conjugate()


def eigen_bounds(H,lo,hi):
    n=H.nrows()
    eye=acb_mat([[int(i==j) for j in range(n)] for i in range(n)])
    certify_positive(H-arb(lo)*eye)
    # 严格负行列式给至少一个小于 hi 的本征值。
    assert (H-arb(hi)*eye).det().real<0


def verify(root):
    lib=run_path(str(root/"ch39-local-yukawa-plucker.py"))
    original=run_path(str(root/"ch39-joint-zero-individual-light-leakage.py"))
    shell=lib["shell"]
    ctx.prec=320
    height=acb.zeta_zero(2_306_888_100_085).imag
    delta=arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    center=15-delta/3
    frames={q:lib["light_basis"](q,height,center) for q in (1,2,-4)}
    for q in frames:
        lap=shell["laplacian"](q,height)
        masses=[abs(z.real*(z.real-center)) for z in lap.eig(multiple=True)]
        assert sum(m<9 for m in masses)==3 and all(m<9 or m>9+(1-delta)**2 for m in masses)
    u=frames[1]
    p=u*u.conjugate().transpose()
    d=[p[i,i].real for i in range(19)]
    assert min(d)>arb(".000026")
    for q in (2,-4):
        v=frames[q]
        target=v*v.conjugate().transpose()
        G=acb_mat([[(int(i==j)-target[i,j])*p[j,i] for j in range(19)] for i in range(19)])
        T=acb_mat([[G[i,j]/(d[i]*d[j]).sqrt() for j in range(19)] for i in range(19)])
        eigen_bounds(T,*CONDITIONAL[q])
        s=arb(".00468" if q==2 else ".00467")
        certify_positive(acb_mat([[G[i,j]-int(i==j)*s*s for j in range(19)] for i in range(19)]))
        # 保迹通道的对偶证书及显式单一对角酉的上界。
        angles=[None]*19
        for group,t in zip(original["GROUPS"],ANGLES[q],strict=True):
            for i in group:angles[i]=arb(t)
        phase=acb_mat([[acb(0,t).exp()] for t in angles])
        output=G*phase
        dual=[(output[i,0]*phase[i,0].conjugate()).real.mid()-EPSILON for i in range(19)]
        certify_positive(acb_mat([[G[i,j]-(dual[i] if i==j else 0) for j in range(19)] for i in range(19)]))
        lo=sum(dual)/3
        hi=(phase.conjugate().transpose()*output)[0,0].real/3
        assert arb(COMPLETE[q][0])<lo<hi<arb(COMPLETE[q][1])
        # 固定一个源态可直接调每站衰减率，将其送到指定的目标轻态。
        assert all(abs(u[i,0])>0 for i in range(19))
        success=arb(1)/10**4
        k=[success.sqrt()*v[i,0]/u[i,0] for i in range(19)]
        assert max(abs(z) for z in k)<1
        print("目标荷、整空间约束满列秩与奇异值下界、后选择条件重率的最优区间",
              q,19,s.str(6),*CONDITIONAL[q],flush=True)
        print("逐站 CPTP 最小平均重率的上下界、单输入无渗漏构造成功率",
              lo.str(15),hi.str(15),success.str(6),flush=True)


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
