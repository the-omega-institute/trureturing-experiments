"""以目标荷原生联络搭建顶点、跨目标荷共享所有复弦距系数的联合重泄漏。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


def adjoint(a):
    return a.conjugate().transpose()


def class_matrices(shell, height, delta, charge):
    pole=shell["soft"]["pole"]
    clocks=shell["geo"]["phases"](height)
    eye=acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    by_distance={}
    for i,j,d,z in shell["geo"]["edges"]:
        spoke=pole in (i,j)
        if spoke:
            amplitude=1-delta if i in shell["soft"]["star"]["orbit"] or j in shell["soft"]["star"]["orbit"] else arb(1)
        else:
            amplitude=(z*clocks.get(d,1))**charge
        key=(int(spoke),d)
        if key not in by_distance:
            by_distance[key]=acb_mat(19,19)
        by_distance[key][i,j]=amplitude
        by_distance[key][j,i]=acb(amplitude).conjugate()
    assert len(by_distance)==23
    return [eye,*[by_distance[k] for k in sorted(by_distance)]]


def gram(blocks):
    a=acb_mat([m.entries() for m in blocks]).transpose()
    return adjoint(a)*a


def verify(root):
    light=run_path(str(root/"ch39-local-yukawa-plucker.py"))
    bounds=run_path(str(root/"ch39-joint-zero-diagonal-light-filters.py"))
    shell=light["shell"]
    ctx.prec=320
    h=acb.zeta_zero(2_306_888_100_085).imag
    delta=arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
    center=15-delta/3
    u=light["light_basis"](1,h,center)
    v=light["light_basis"](2,h,center)
    w=light["light_basis"](-4,h,center)
    overlap=adjoint(v)*w
    Qs=[f*adjoint(f) for f in (v,w)]
    heavy=acb_mat(24,24)
    total=acb_mat(24,24)
    separate=[]
    for q,Q in zip((2,-4),Qs,strict=True):
        basis=class_matrices(shell,h,delta,q)
        out=[m*u for m in basis]
        leak=[z-Q*z for z in out]
        H=gram(leak)
        D=gram(out)
        heavy+=H
        total+=D
        separate.append((H,D))
    # 相同输出时，对两目标均匀取平均的泄漏至少为 (1-σmax)/2。
    principal=adjoint(overlap)*overlap
    eig=[t.real for t in principal.eig(algorithm="rump")]
    largest=max(eig)
    assert arb(".99870222")<largest.sqrt()<arb(".99870224")
    principal_floor=(1-largest.sqrt())/2
    assert principal_floor>arb(".00064888")
    print("两目标轻空间的最大主角余弦与同输出平均泄漏锐下界",
          largest.sqrt().str(14),principal_floor.str(14),flush=True)
    positive=bounds["certify_positive"]
    identity=acb_mat([[int(i==j) for j in range(24)] for i in range(24)])
    positive(heavy-arb("8.5e-8")*identity)
    lo,hi=arb(".00015015280"),arb(".00015015282")
    positive(heavy-lo*total)
    assert (heavy-hi*total).det().real<0
    print("共享复系数重输出秩与联合条件重泄漏锐下界区间",24,lo.str(12),hi.str(12),flush=True)

    # 对两目标分别条件化后取较差者。加权 Gram 铅笔给下界，其最小本征方向给上界。
    t=arb(".49301798313")
    weighted_H=t*separate[0][0]+(1-t)*separate[1][0]
    weighted_D=t*separate[0][1]+(1-t)*separate[1][1]
    low=arb(".00015018063")
    estimate=arb(".000150180638")
    positive(weighted_H-low*weighted_D)
    zl=weighted_H-estimate*weighted_D
    indices=[i for i in range(24) if i!=22]
    reduced=acb_mat([[zl[i,j] for j in indices] for i in indices])
    assert not reduced.det().contains(0)
    solution=reduced.solve(acb_mat([[-zl[i,22]] for i in indices])).entries()
    solution.insert(22,acb(1))
    trial=acb_mat([[c] for c in solution])
    achieved=[]
    for H,D in separate:
        heavy_part=(adjoint(trial)*H*trial)[0,0].real
        output=(adjoint(trial)*D*trial)[0,0].real
        assert output>0
        achieved.append(heavy_part/output)
    cap=arb(".00015018065")
    assert all(rate<cap for rate in achieved),achieved
    print("两目标中较差条件泄漏的最优严格区间",low.str(14),cap.str(12),flush=True)


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
