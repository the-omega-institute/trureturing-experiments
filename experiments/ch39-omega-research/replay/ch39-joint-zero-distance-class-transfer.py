"""原生荷一联络按无向弦距共系数的整代传递及弦距多项式阶数阈值。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


def adjoint(a):
    return a.conjugate().transpose()


def columns(outputs):
    return acb_mat([a.entries() for a in outputs]).transpose()


def shell_classes(root):
    light=run_path(str(root/"ch39-local-yukawa-plucker.py"))
    cyclic=run_path(str(root/"ch39-joint-zero-cyclic-light-transfer.py"))
    shell=light["shell"]
    ctx.prec=320
    h=acb.zeta_zero(2_306_888_100_085).imag
    delta=arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
    S,edges=cyclic["generator"](shell)
    center=15-delta/3
    eye=acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    pole=shell["soft"]["pole"]
    phases=shell["geo"]["phases"](h)
    outer,spokes={},{}
    for i,j,length,z in shell["geo"]["edges"]:
        if pole in (i,j):
            store=spokes
            w=1-delta if i in shell["soft"]["star"]["orbit"] or j in shell["soft"]["star"]["orbit"] else arb(1)
        else:
            store=outer
            w=z*phases.get(length,1)
        if length not in store:
            store[length]=acb_mat(19,19)
        store[length][i,j]=w
        store[length][j,i]=acb(w).conjugate()
    ds=sorted(outer); refs=[eye]+[outer[d] for d in ds]
    assert len(ds)==20 and sorted(spokes)==[68,136,544]
    allclasses=[*refs,*spokes.values()]
    assert all(all(z.contains(0) for z in (S*A-A*S).entries()) for A in allclasses)
    u=light["light_basis"](1,h,center)
    order=cyclic["character_order"](u,S)
    assert sorted(order)==[0,1,2]
    return light, cyclic, S, h, center, eye, ds, refs, spokes, u, order


def verify(root):
    light, cyclic, S, h, center, eye, ds, refs, spokes, u, order = shell_classes(root)
    outer = dict(zip(ds, refs[1:], strict=True))
    positive=run_path(str(root/"ch39-joint-zero-diagonal-light-filters.py"))["certify_positive"]
    # 选择十九根独立的受限列。十九维是全部从正则 C3 源表示到十九站的交织空间维数。
    indices=(15,16,10,6,9,14,18,0,8,4,1,17,5,7,13,11,12,20,3)
    select=[refs[i] for i in indices]
    M=columns([a*u for a in select])
    gram=adjoint(M)*M
    positive(gram)
    polynomials=[eye]+[spokes[k] for k in sorted(spokes)]
    # 只在非辐边上要求一个多项式。横坐标是整弦长平方 d。
    ts={d:[arb(1),arb(d)/400-1] for d in ds}
    for d,row in ts.items():
        while len(row)<13:
            row.append(2*row[1]*row[-1]-row[-2])
    for n in range(13):
        polynomials.append(sum((ts[d][n]*outer[d] for d in ds),acb_mat(19,19)))
    A=columns([p*u for p in polynomials[:16]])
    total=adjoint(A)*A
    for q,limits in ((2,("3.96e-9","3.97e-9")),(-4,(".0001494",".0001495"))):
        v0=light["light_basis"](q,h,center)
        vc=cyclic["character_order"](v0,S)
        v=acb_mat([[v0[i,vc.index(k)] for k in order] for i in range(19)])
        Q=eye-v*adjoint(v)
        coeff=gram.solve(adjoint(M)*columns([v]))
        T=sum((k*a for k,a in zip(coeff.entries(),select,strict=True)),acb_mat(19,19))
        assert all(z.contains(0) and abs(z)<arb("1e-32") for z in (T*u-v).entries())
        norm=sum(abs(z)**2 for z in T.entries()).sqrt()
        assert norm<1000
        print("荷、全部弦距模型的输出秩与重输出秩、等距拟合的 Frobenius 范数、整代统一成功率",
              q,19,16,norm.str(13),".000001",flush=True)
        heavy=[Q*a*u for a in polynomials]
        G=adjoint(columns(heavy[:16]))*columns(heavy[:16])
        lo,hi=map(arb,limits)
        positive(G-lo*total)
        assert (G-hi*total).det().real<0
        # 第十二次项的重输出位于同一十六维交织空间，首十六项恰给一组基。
        last=columns([heavy[16]])
        fit=G.solve(adjoint(columns(heavy[:16]))*last)
        T=polynomials[16]-sum((k*a for k,a in zip(fit.entries(),polynomials[:16],strict=True)),acb_mat(19,19))
        out=T*u
        errors=(Q*out).entries()
        assert all(z.contains(0) and abs(z)<arb("1e-6") for z in errors), max(abs(z).upper() for z in errors)
        restriction=adjoint(v)*out
        assert min(abs(restriction[i,i]) for i in range(3))>0
        assert all(restriction[i,j].contains(0) for i in range(3) for j in range(3) if i!=j)
        scale=sum(abs(z)**2 for z in T.entries()).sqrt().upper()
        success=min(abs(restriction[i,i]).lower()**2 for i in range(3))/scale**2
        assert success>arb("1e-12")
        print("十一次的最优条件泄漏区间、十二次零漏三线传输的统一成功率下界",
              *limits,success.lower().str(10),flush=True)


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
