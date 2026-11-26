"""实际零点的两步局部质量块与全部中间阈值的两圈规范交点。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


ROOT=Path(__file__).parent


def active_phases(sm,weights,light):
    q=sm["q"]
    b=sm["BETA"]
    # 每个活动模有三份 Dirac。|q|=3 时另有一个复标量弱双重态。
    B=sm["B"]
    beta=b
    phases=[(beta,B)]
    for _,c in light:
        color,weak=(1,2) if c==3 else (3,1) if c==4 else (3,2)
        index=[arb(color*weak*c*c)/60,arb(color)/2 if weak==2 else 0,arb(weak)/2 if color==3 else 0]
        casimir=[arb(c*c)/60,q("3/4") if weak==2 else 0,q("4/3") if color==3 else 0]
        scalar=int(c==3)
        add=[[index[i]*((12+4*scalar)*casimir[j]+
              ((20+arb(2)*scalar/3)*(0,2,3)[i] if i==j else 0))
              for j in range(3)] for i in range(3)]
        beta=[v+k for v,k in zip(beta,weights[c])]
        B=[[v+k for v,k in zip(r,dr)] for r,dr in zip(B,add)]
        phases.append((beta,B))
    return phases


def verify():
    model,_,rows,initial,weights=run_path(str(ROOT/"ch39-neutral-mean-local-mass.py"))["inputs"]()
    running=run_path(str(ROOT/"ch39-two-tier-dirac-running.py"))["integrate_segments"]
    sm=run_path(str(ROOT/"ch39-two-loop-shell-separation.py"))
    ctx.prec=192
    light=rows[:6]
    assert [c for _,c in light]==[1,1,3,1,4,3]
    assert all(x[0]<y[0] for x,y in zip(light,light[1:]))
    b=sm["BETA"]
    matrix=arb_mat([[b[i],sum(weights[c][i] for _,c in light),2*arb.pi()] for i in range(3)])
    base=[2*arb.pi()*a+sum(weights[c][i]*v.log() for v,c in light) for i,a in enumerate(initial)]
    phases=active_phases(sm,weights,light)

    def times_image(L,t,steps):
        assert (rows[5][0].log()<t<rows[6][0].log())
        times=[L-t+light[0][0].log()]
        times.extend((w[0]/v[0]).log() for v,w in zip(light,light[1:]))
        times.append(t-light[-1][0].log())
        assert all(v>0 for v in times)
        segments=[(length,bi,Bi,steps if k==0 else 24)
                  for k,(length,(bi,Bi)) in enumerate(zip(times,phases))]
        corr,_=running(initial,segments,sm["D"],sm["C"])
        return matrix.solve(arb_mat([[v+2*arb.pi()*w] for v,w in zip(base,corr)])).entries()
    L,t=arb("29.57").union(arb("29.96")),arb("-1.21").union(arb("-1.00"))
    middle=times_image(L,t,512)
    assert L.contains_interior(middle[0]) and t.contains_interior(middle[1])
    for _ in range(4):
        L,t=L.intersection(middle[0]),t.intersection(middle[1])
        middle=times_image(L,t,2048)
    L,t=L.intersection(middle[0]),t.intersection(middle[1])
    common=middle[2]
    mu=model["mz"]*L.exp()
    unit=model["mz"]*(L-t).exp()
    assert arb("7.5e14")<mu<arb("7.9e14") and arb("2.20e15")<unit<arb("2.45e15")
    assert 38<common<39
    print("两圈尖阈值交点、质量单位、共同逆耦合的包围",
          *[(v.lower().str(12),v.upper().str(12)) for v in (mu,unit,common)])
    print("六个阈值的质量区间 / GeV",
          [(float((v*unit).lower()),float((v*unit).upper())) for v,_ in light])


if __name__=="__main__":
    verify()
