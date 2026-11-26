"""三枚实际零点的共同截止，同时规定中间与终点的正中性动能。"""
from pathlib import Path
from runpy import run_path
from flint import arb,arb_mat,ctx


def verify(root):
    gravity=run_path(str(root/"ch39-three-index-free-atoms.py"))
    moment,_=gravity["spectra"](root)
    gauge=run_path(str(root/"ch39-all-field-signed-crossing.py"))
    neutral=run_path(str(root/"ch39-full-neutral-free-cutoff.py"))
    rotation=run_path(str(root/"ch39-rotation-portal-top-energy.py"))
    signs,masses=rotation["rotation_signs"](root,2306888100085)
    data=gauge["joint_system"](root,"1.176397e14",order=12,last=36)
    ctx.prec=320
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    portal=arb(".42995401872")
    species,light=neutral["response_rows"](data,masses,signs,rho,portal)
    sites=list(map(arb,("1e-50","4e-15","1e-12","3.5e-11","1e-7","1.23e-6",".5","1",
                        "1.7033e-9","2.54e-7","8.96948e-6",".0303688",
                        "1.72782598e-20","2.44061907e-10","7.28618175e-17","8.65964323e-12")))
    rule=gauge["angular_rule"](18,60)
    p0=data["p0"];p1=data["p1"];pstar=(arb("1e12")/arb("1.22089e19"))**2
    pi=arb.pi();heat=data["heat"]
    def gauge_rows(x,p):
        result=[arb(0)]*3
        log=(p/p0).log()
        for kind,z,beta in species:
            term=gauge["angular_exposure"](kind,z,x,p0,p,log,rule)/(4*pi)
            for j,b in enumerate(beta):
                result[j]-=b*term
        return result
    def potential(x):
        slope=sum(g*k*x*heat(2,k*rho/x)/4 for g,k in light)
        curv=-sum(g*k**2*heat(1,k*rho/x)/4 for g,k in light)
        for z,s in zip(masses,signs):
            arg=(z+s*portal*rho)/x
            slope+=s*portal*x*heat(2,arg)
            curv-=portal**2*heat(1,arg)
        return [slope/(8*pi**2),curv/(8*pi**2)]
    cols=[]
    for x in sites:
        qq,qt,colour=gauge_rows(x,p1)
        cols.append(moment(x)+potential(x)+[qq-8*colour/3,qt-colour,neutral["weak_column"](x,p0,p1),
                    colour,*gauge_rows(x,pstar)])
    a1,a2,a3=data["initial"]
    em=5*a1/3+a2
    target=arb_mat([[v] for v in (1,1,arb("1e-122"),0,0,0,0,0,arb(".258"),8*a3/3-em,a3-a2,
                                    0,39-a3,104-em,39-a2,39-a3)])
    matrix=arb_mat(cols).transpose()
    mid=arb_mat([[v.mid() for v in row] for row in matrix.tolist()])
    estimate=mid.solve(target)
    centers=arb_mat([[v.mid()] for v in estimate.entries()])
    pre=mid.inv()
    uncertainty=arb_mat([[abs(v).upper() for v in row] for row in (pre*matrix-arb_mat([[int(i==j) for j in range(16)] for i in range(16)])).tolist()])
    residual=pre*(target-matrix*centers)
    widths=[abs(v).upper()+arb("1e-70") for v in residual.entries()]
    for _ in range(12):
        defects=(uncertainty*arb_mat([[u] for u in widths])).entries()
        widths=[(abs(v).upper()+e)*arb("1.01")+arb("1e-70") for v,e in zip(residual.entries(),defects)]
    end=(uncertainty*arb_mat([[u] for u in widths])).entries()
    assert all(abs(v)+e<rad for v,e,rad in zip(residual.entries(),end,widths))
    weights=[v+arb(0,rad.upper()) for v,rad in zip(centers.entries(),widths)]

    assert all(not v.contains(0) for v in weights)
    running=run_path(str(root/"ch39-full-neutral-running.py"))
    weak=running["weak_bound"](data,sites,weights,run_path(str(root/"ch39-neutral-axial-top.py")),
                               run_path(str(root/"ch39-neutral-bosonic.py")))
    cost=sum(-v for v in weights if v<0)
    print("十六权负量、全程弱差绝对上界",cost.str(13),weak.str(10),flush=True)
    print("十六权",*[v.str(10) for v in weights],flush=True)
    assert weak<arb("4e-11") and arb("3.407e7")<cost<arb("3.410e7")
    # 低能段先提取与截止无关的全物种核；大质量行可直接按 1/z 估计。
    p_low=(arb("7e7")/arb("1.22089e19"))**2
    total=sum(abs(c).upper() for c in weights)
    invx=sum(abs(c).upper()/x for x,c in zip(sites[1:],weights[1:]))
    drift=[arb(0)]*3
    for kind,z,beta in species:
        if z<arb("1e-25"):
            bound=abs(1-weights[0])*(p_low/p0).log()+invx*(p_low-p0)/4
            bound+=abs(weights[0])*12*sites[0]/p0
        else:
            bound=total*(p_low-p0)/(4*z)
        for j,coef in enumerate((beta[1],beta[0]-beta[1],beta[2])):
            drift[j]+=abs(arb(coef)).upper()*bound/(4*pi)
    lows=[a-b for a,b in zip((a2,em-a2,a3),drift)]
    print("从 M_Z 到 7e7 GeV 的 QT3、QQ−QT3、色下界",*[a.lower().str(14) for a in lows],flush=True)
    assert min(v.lower() for v in lows)>1 and lows[0]*lows[1]/(lows[0]+lows[1])>weak
    # 反例读数与终点作差并保留精确归一，避免百万量级权重的不相关重复误差。
    at_bad=(arb("1e13")/arb("1.22089e19"))**2
    bad=[gauge_rows(x,at_bad)[0] for x in sites]
    bad_shift=[bad[i]-(cols[i][9]+8*cols[i][12]/3) for i in range(len(sites))]
    ref=bad_shift[7].mid()
    qq_bad=104+ref+sum(c*(v-ref) for c,v in zip(weights,bad_shift))
    print("1e13 GeV 的电磁逆耦合",qq_bad.str(12),flush=True)
    assert arb("-2.83e5")<qq_bad<arb("-2.76e5")
    # 关于 log Q 的一阶导数，以正角核归一和 h p/(z+hp)≤1 直接控制。
    qSlope=[sum(abs(arb(b[j])).upper() for kind,z,b in species if b[j]!=0)
            *sum(abs(v).upper() for v in weights)/(2*pi) for j in range(3)]
    common=max(qSlope[0].upper()+qSlope[1].upper(),qSlope[1].upper(),qSlope[2].upper())
    h=arb("2e-9")
    print("两个能标附近的 log Q 半宽与三行最大漂移",h.str(4),(common*h).upper().str(12),flush=True)
    b,y,colour=arb(29),arb(55),arb(29)
    assert colour>0 and b*y/(b+y)>weak
    assert common*h<10
    return data,sites,weights,species


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
