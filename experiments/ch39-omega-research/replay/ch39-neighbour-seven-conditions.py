"""同一七条件解的相邻实际零点响应，以及两点共同匹配的必要条件。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def dot(weights, row):
    return sum(w*d for w, d in zip(weights,row,strict=True))


def verify():
    root = Path(__file__).resolve().parent
    build = run_path(str(root / "ch39-six-broken-conditions.py"))["broken_system"]
    ctx.prec = 192
    energy = arb("1.1763969e14").union(arb("1.17639745e14"))
    portal = arb(".42995401871").union(arb(".42995401873"))
    n = 2_306_888_100_085
    ref = build(root, energy, portal, number=n)
    c, _, _ = ref["weights_at"](portal)
    mat0, _ = ref["columns_at"](portal)

    changes = []
    for other in (n+1,n-1):
        adj = build(root, energy, portal, number=other)
        mat1, _ = adj["columns_at"](portal)
        diff = mat1-mat0
        shifts = [dot(c,[diff[i,j] for j in range(5)]) for i in range(1,5)]
        # 三循环两种符号的总数相同。先逐节点消去常数 1，保留微小的谱效应。
        heat = ref["original"]["heat"]
        assert sum(ref["signs"]) == sum(adj["signs"])
        def heavy_slope(model, x):
            terms = []
            for sign,z in zip(model["signs"],model["masses"],strict=True):
                mass = z+sign*portal*model["rho"]
                terms.append(sign*(x*(-mass/x).expm1()-mass*heat(1,mass/x)))
            return sum(terms)
        slope = portal*sum(w*(heavy_slope(adj,x)-heavy_slope(ref,x))
                           for w,x in zip(c,ref["original"]["sites"],strict=True))/(8*arb.pi()**2)
        sym0, sym1 = ref["original"], adj["original"]
        mixed = dot(c,[(v[4]-v[5])-(u[4]-u[5]) for u,v in zip(sym0["cols"],sym1["cols"],strict=True)])
        # 两端的轻 W 移质完全相同。各自重标量的带号移质用导数的一致界控制。
        scalar_error = (ref["rho"]*sum(abs(w) for w in c)/(24*arb.pi()) *
                        sum(1/(z-ref["rho"]) for d in (ref,adj) for z in d["masses"]))
        mixed += arb(0,scalar_error.upper())
        print(other, "第 K,W,B,曲率,斜率,C 的位移",*[v.str(12) for v in (*shifts,slope,mixed)], flush=True)

        changes.append(shifts[:2])
        kwb, vwb, bfit, curv = shifts
        if other > n:
            assert arb("1.05657e-7")<kwb<arb("1.05658e-7") and arb("2.27663e-8")<vwb<arb("2.27664e-8")
            assert arb(".0038")<bfit<arb(".0042") and arb(".0011384")<curv<arb(".0011386")
            assert arb(".0024")<mixed<arb(".0036")
            assert -arb("1.438e-12") < slope < -arb("1.437e-12")
        else:
            assert -arb("6.07544e-8")<kwb<-arb("6.07542e-8") and -arb("1.43663e-8")<vwb<-arb("1.43662e-8")
            assert -arb(".0056")<bfit<-arb(".0044") and -arb(".0010570")<curv<-arb(".0010568")
            assert -arb(".0030")<mixed<-arb(".0028")
            assert arb("3.376e-12") < slope < arb("3.378e-12")

    # 给邻点总概率一时，最优比例令 K、W 偏差互为相反数。
    (kp,wp), (km,wm) = changes
    frac = -(km+wm)/(kp-km+wp-wm)
    minimum = abs((km*wp-kp*wm)/(kp-km+wp-wm))
    assert 0 < frac < 1 and kp>0>km and wp>0>wm
    assert arb(".36906") < frac < arb(".36907")
    assert arb("6.6196e-10") < minimum < arb("6.6200e-10")
    print("相邻编号全部质量分给两端时，右端最优份额、引力最大偏差的最小值",frac.str(12),minimum.str(12))


if __name__ == "__main__":
    verify()
