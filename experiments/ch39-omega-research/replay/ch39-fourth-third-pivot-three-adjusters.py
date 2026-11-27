"""固定第四壳第三枢轴，以三个低截止补偿其他十六个高截止的任意取值。"""
from pathlib import Path
from runpy import run_path
from math import factorial

from flint import arb, arb_mat, ctx


def local_box():
    """按第 10、1、17 层返回截止中心与对数半宽，高层按原层号递增。"""
    mid=list(map(arb,(
        ".0004565202895202273",".008089448084270673",".013037841943892032",".0006540513601246211",
        ".00569127680411594",".0007318857601587009",".00012052475355810329",".0016832907773903",
        ".0008412347835513889","9.33255615265773e-5",".001430388180319893","9.450157091069868e-5",
        ".09084640056401722",".019211171771208315",".01694232219822086",".015156506382214815")))
    return dict(layers=(9,0,16), centers=list(map(arb,("4.156672262801198e-34",
                "9.079019273084664e-10","1.355984621446739e-23"))),
                radii=list(map(arb,("9e-10","4e-9","1e-8"))),
                high_centers=mid, high_log_radius=arb(".01"))


def verify(root, *, number=2_306_888_100_089, box_data=None):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    native=lib("native-rh-ray-pvm")
    _,_,pivots,_=next(row for row in native["shell_weights"](root,number) if row[0]==6)
    c=[a.real for a in pivots[2]]
    gauge=lib("all-field-signed-crossing")
    data=gauge["joint_system"](root,"1.176397e14",number=number,order=4,last=8)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec=192
    neutral=lib("full-neutral-free-cutoff")
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    species,_=neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    pi=arb.pi()
    rows=[(k,z,[(qq-8*col/arb(3))/(4*pi),(mix-col)/(4*pi)]) for k,z,(qq,mix,col) in species]
    a1,a2,a3=data["initial"]
    starts=[5*a1/3+a2-8*a3/3,a2-a3]
    p0,p1=data["p0"],data["p1"]
    rule=gauge["angular_rule"](20,55)
    log=(p1/p0).log()
    def exposures(x):
        return [gauge["angular_exposure"](k,z,x,p0,p1,log,rule) for k,z,_ in rows]
    floor=arb("1e-6")
    ref=arb(".5")
    hmin,href=exposures(floor),exposures(ref)
    order=22
    beta=[arb(factorial(i))**2/factorial(2*i+1) for i in range(order+2)]
    angular=((0,6,1),(3,-12,1),(5,12,7),(9,12,11))
    polys=[[(a*beta[j]+b*beta[j+1])*(-1)**j/(d*factorial(j)) for j in range(order+1)] for a,b,d in angular]
    def avg(t,pol):
        if t>100:
            return arb(0).union(12/t)
        assert t<1
        v=arb(0)
        for a in reversed(pol):
            v=v*t+a
        return v+arb(0,((t/4)**(order+1)/factorial(order+1)).upper())
    def avgband(t,pol):
        return avg(t.upper(),pol).union(avg(t.lower(),pol))
    def deriv(x):
        fac=[avgband(p0/x,pol)-avgband(p1/x,pol) for pol in polys]
        return [-sum(b[j]*(-z/x).exp()*fac[k] for k,z,b in rows) for j in range(2)]
    weak=native["weak_series"](56)
    local=local_box() if box_data is None else box_data
    tcenter=local["centers"][0]
    high_mid=local["high_centers"]
    width=local["radii"][0]
    last=16
    modes = [
            ("local",*local["centers"][1:],*local["radii"][1:]),
            ("9e-5","9.079019140623481e-10","1.355984336417077e-23",".00035",".000075"),
            ("1e-6","9.079019140623481e-10","1.355984336417077e-23",".035",".0065")]
    for mode,m1,m2,r1,r2 in modes if box_data is None else modes[:1]:
        active=(0,last)
        other=[i for i in range(19) if i not in (9,*active)]
        center=[arb(m1),arb(m2)]
        base=[exposures(x) for x in (tcenter,*center)]
        refs=[exposures(x) for x in high_mid] if mode=="local" else [href]*16
        const=[starts[j]-sum(b[j]*(sum(c[i]*e[k] for i,e in zip((9,*active),base,strict=True))
                                +sum(c[i]*e[k] for i,e in zip(other,refs,strict=True))) for k,(_,_,b) in enumerate(rows))
               for j in range(2)]
        def jac(v):
            return arb_mat([[c[i]*a for a in deriv(x*u.exp())] for i,x,u in zip(active,center,v,strict=True)]).transpose()
        pre=jac([arb(0)]*2).mid().inv().mid()
        assert not pre.det().contains(0)
        # 高端单核单调；从参照点到无穷大的差不超过 (p1-p0)/(4 ref)。
        # 将每个系数向量先左乘 pre 抑制分量间的依赖，再按高层符号求和。
        projected=[[sum(pre[j,k]*b[k] for k in range(2)) for _,_,b in rows] for j in range(2)]
        off=pre*arb_mat([[v] for v in const])
        if mode=="local":
            high_width=local["high_log_radius"]
            high_jac=[deriv(x*arb(0,high_width.upper()).exp()) for x in high_mid]
        else:
            hmin=exposures(arb(mode))
        tgrad=pre*arb_mat([[a*c[9]] for a in deriv(tcenter*arb(0,width.upper()).exp())])
        for j,b in enumerate(projected):
            exposure_delta=[(l-r).union((p1-p0)/(4*ref)) for l,r in zip(hmin,href,strict=True)]
            perturb=-sum(v*delta for v,delta in zip(b,exposure_delta,strict=True))
            if mode=="local":
                high_bound=sum(abs(c[i]*sum(pre[j,k]*g[k] for k in range(2))) for i,g in zip(other,high_jac,strict=True))
                off[j,0]+=arb(0,(high_width*high_bound).upper())
            else:
                off[j,0]+=sum(c[i]*perturb for i in other)
            off[j,0]+=arb(0,(abs(tgrad[j,0])*width).upper())
        radii=[arb(r1),arb(r2)]
        box=[arb(0,r.upper()) for r in radii]
        defect=arb_mat([[1,0],[0,1]])-pre*jac(box)
        ratios=[(abs(off[j,0])+sum(abs(defect[j,k])*radii[k] for k in range(2)))/radii[j] for j in range(2)]
        if mode=="local":
            local["preconditioner"]=pre
        print(mode,"高端一致的预条件位移、盒半径、压缩自映上界",last+1,[v.str(12) for v in off.entries()],
              [v.str(5) for v in radii],[v.upper().str(12) for v in ratios],flush=True)
        assert all(v<1 for v in ratios)
        lows=[v*u.exp() for v,u in zip(center,box,strict=True)]
        # 高截止弱差包围半径随 x 单调减小。在两个参数端点分别处理唯一极低层。
        highweak=neutral["weak_column"](floor,p0,p1)
        rest=sum(c[i]*highweak for i in other)+sum(c[i]*neutral["weak_column"](x,p0,p1)
                            for i,x in zip(active,lows,strict=True))
        ends=[rest+c[9]*weak(tcenter*(s*width).exp()/p0) for s in (-1,1)]
        assert ends[0]>0>ends[1]
        print("弱差两端与三层截止盒",
              [(v.lower().str(12),v.upper().str(12)) for v in ends],
              [(v.lower().str(15),v.upper().str(15)) for v in [tcenter*arb(0,width.upper()).exp(),*lows]],flush=True)

    return local


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
