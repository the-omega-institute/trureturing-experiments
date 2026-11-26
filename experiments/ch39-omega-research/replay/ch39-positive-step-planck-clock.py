"""正的单步传递与二进分母，固定普朗克归一后从两圈交点反求强耦合。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat, ctx


ROOT=Path(__file__).parent


def positive_denominators(phase):
    cut=phase["cut"]
    signless=acb_mat(19,19)
    for i,j,_,_ in cut["SHELL"]["edges"]:
        signless[i,i]+=1
        signless[j,j]+=1
        signless[i,j]=signless[j,i]=1
    assert cut["count_above"](signless,arb("30.38578"))==1
    assert cut["count_above"](signless,arb("30.38580"))==0
    for q in (1,2,3,4,6):
        assert phase["count_below"](q,arb("22.9999"))==19
    assert phase["count_below"](1,arb("22.1"))<19
    print("任意单位边相位、原生钟相位、共同的最小二次幂分母",31,23,32)


def finite_match(model,rows,weights,initial,den,sharp):
    pi=arb.pi()
    unit=arb("2.435e18")/den**2
    def intercepts(L):
        out=[2*pi*v-b*L for v,b in zip(initial,model["beta"])]
        for h,q in rows:
            mass=(unit*h/model["mz"])**2
            def exposure(kernel):
                rlo=mass.lower()*(-2*L.upper()).exp()
                rhi=mass.upper()*(-2*L.lower()).exp()
                val=kernel(rhi,arb(1)).union(kernel(rlo,arb(1)))
                return val-arb(0).union(1/(8*mass.lower()))
            if not q:
                v=exposure(lambda x,y:model["sheet"](x,y,True))
                out=[a-k*v for a,k in zip(out,weights[0])]
            else:
                f=exposure(lambda x,y:model["sheet"](x,y,False))
                out=[a-6*k*f for a,k in zip(out,model["weyl"][model["charges"].index(q)])]
                if q==3:
                    hval=exposure(model["scalar"])
                    out=[a-k*hval for a,k in zip(out,model["higgs"])]
        return out

    y,w,c=intercepts(sharp[0])
    print("在尖阈值一圈原交点上采用有限核后的 Y-W、C-W 逆耦合差",den,
          ((y-w)/(2*pi)).str(10),((c-w)/(2*pi)+sharp[2]-initial[2]).str(10))
    lower,upper=arb("31.60"),arb("31.88")
    def difference(t):
        a,b,_=intercepts(t)
        return a-b
    assert difference(lower)<0<difference(upper)
    for _ in range(21):
        mid=(lower+upper)/2
        d=difference(mid)
        if d>0:
            upper=mid
        else:
            assert d<0
            lower=mid
    L=lower.union(upper)
    y,_,c=intercepts(L)
    alpha=1/(initial[2]+(y-c)/(2*pi))
    assert arb(".1531")<alpha<arb(".158")
    print("有限动量一圈另一匹配处方的中间交点",den,
          (model["mz"]*L.exp()).str(9),alpha.str(9))


def verify():
    phase=run_path(str(ROOT/"ch39-neutral-mean-phase-inertia.py"))
    ctx.prec=192
    positive_denominators(phase)
    model,_,rows,initial,weights=run_path(str(ROOT/"ch39-neutral-mean-local-mass.py"))["inputs"]()
    flow=run_path(str(ROOT/"ch39-two-tier-dirac-running.py"))["integrate_segments"]
    sm=run_path(str(ROOT/"ch39-two-loop-shell-separation.py"))
    phases=run_path(str(ROOT/"ch39-neutral-mean-two-loop.py"))["active_phases"]
    ctx.prec=192
    light=rows[:6]
    assert [c for _,c in light]==[1,1,3,1,4,3]
    regimes=phases(sm,weights,light)
    pi=arb.pi()
    for den,seed in [(31,".12103"),(32,".11889")]:
        kappa=arb("2.435e18")/den**2
        knots=[(kappa*v/model["mz"]).log() for v,_ in light]
        next_knot=(kappa*rows[6][0]/model["mz"]).log()
        total=regimes[-1][0]
        matrix=arb_mat([[total[i],2*pi,-2*pi*int(i==2)] for i in range(3)])
        base=[(2*pi*initial[i] if i!=2 else 0)+sum(weights[q][i]*r for (_,q),r in zip(light,knots))
              for i in range(3)]
        sharp=matrix.solve(arb_mat([[a] for a in base])).entries()
        assert knots[-1]<sharp[0]<next_knot
        def image(L,inverse,steps):
            assert knots[-1]<L<next_knot
            durations=[knots[0],*[b-a for a,b in zip(knots,knots[1:])],L-knots[-1]]
            segments=[(t,b,B,steps if k==0 else max(48,steps//32))
                      for k,(t,(b,B)) in enumerate(zip(durations,regimes))]
            corr,_=flow([*initial[:2],inverse],segments,sm["D"],sm["C"])
            return matrix.solve(arb_mat([[a+2*pi*c] for a,c in zip(base,corr)])).entries()
        interval_L=arb("29.70").union(arb("30.03"))
        inverse=1/arb(seed)+arb(0,".12")
        result=image(interval_L,inverse,1024)
        assert interval_L.contains_interior(result[0]) and inverse.contains_interior(result[2])
        for steps in (1024,2048,4096,8192,16384):
            interval_L=interval_L.intersection(result[0])
            inverse=inverse.intersection(result[2])
            result=image(interval_L,inverse,steps)
        L=interval_L.intersection(result[0])
        alpha=1/inverse.intersection(result[2])
        mu=model["mz"]*L.exp()
        print("分母与中间电弱交点的一圈强耦合",den,(1/sharp[2]).str(12))
        print("两圈交点 / GeV 与反推强耦合",*[(v.lower().str(11),v.upper().str(11)) for v in (mu,alpha)])
        finite_match(model,rows,weights,initial,den,sharp)
        if den==32:
            assert arb(".1188")<alpha<arb(".119")
        else:
            assert arb(".1209")<alpha<arb(".1212")


if __name__=="__main__":
    verify()
