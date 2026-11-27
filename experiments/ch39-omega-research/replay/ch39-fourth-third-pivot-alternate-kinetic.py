"""第四壳第三枢轴的另一组五式解及其独立能标动能检验。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    n=2_306_888_100_085
    _,ordered=lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=n+1)[2]
    source=[v.sqrt() for v in ordered]
    native=lib("native-rh-ray-pvm")
    _,_,pivots,_=next(row for row in native["shell_weights"](root,n+4) if row[0]==6)
    c=[v.real for v in pivots[2]]
    assert sum(c).contains(1)
    moments,heat,(field,)=lib("three-index-free-atoms")["broken_fields"](root,[n+4])
    gauge=lib("all-field-signed-crossing")
    mz=arb("91.1876")
    energy=arb("1.176397e14")
    data=gauge["joint_system"](root,energy,number=n+4,order=4,last=8)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,n+4)
    ctx.prec=192
    p0,p1=data["p0"],data["p1"]
    assert p1.overlaps(p0*(energy/mz)**2)
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    species,_=lib("full-neutral-free-cutoff")["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    rows=[(kind,z,(qq-arb(8)*colour/3,mix-colour,colour)) for kind,z,(qq,mix,colour) in species]
    a1,a2,a3=data["initial"]
    starts=(5*a1/3+a2-8*a3/3,a2-a3)
    pi=arb.pi()


    # x 为截止平方。固定十三层，第十层作为外部弱差参数，其余五个自由层按乘法坐标求解。
    anchors=list(map(arb,(
        "8.32787465686893e-11", "1.29676588455162e-17", "6.03963089553863e-25", "0.116010852860954",
        "5.43492569418729e-17", "4.50521660063842e-28", "6.91887563115565e-15", "0.00222337855556152",
        "0.161480494470546", "3.30777788299359e-34", "0.000300447229307098", "0.000132294576646989",
        "8.7162619133473e-11", "0.057757860092149", "0.0013178606810481", "6.15268844055002e-24",
        "0.00810011348202119", "2.63121303105515e-22", "2.02261842030907e-33")))
    selected=(17, 16, 3, 8, 0)
    center=[anchors[i] for i in selected]
    width=arb('3e-8')
    rule=gauge["angular_rule"](19,54)
    log=(p1/p0).log()
    def values(x):
        col=[arb(0)]*3
        for kind,z,beta in rows:
            exposure=gauge["angular_exposure"](kind,z,x,p0,p1,log,rule)/(4*pi)
            for j in range(3):
                col[j]-=beta[j]*exposure
        return [*col[:2],*moments(field,heat,x),col[2]]
    # 对 log x 的规范导数使用指数的 beta 矩。t=p1/x 在五个变量的小区间上有界。
    order=80
    beta=[arb(factorial(j))**2/factorial(2*j+1) for j in range(order+2)]
    angular=((0,6,1),(3,-12,1),(5,12,7),(9,12,11))
    polys=[[(a*beta[j]+b*beta[j+1])*(-1)**j/(d*factorial(j)) for j in range(order+1)] for a,b,d in angular]
    def factors(x):
        def avg(t,poly):
            # 对称半区 h=u(1-u)>=u/2，任一归一化角权不超过 3。
            # 故两个半区的指数积分合计至多 12/t。
            if t>100:
                return arb(0).union(12/t)
            v=arb(0)
            for coeff in reversed(poly):
                v=v*t+coeff
            return v+arb(0,((t/4)**(order+1)/factorial(order+1)).upper())
        return [avg(p0/x,poly)-avg(p1/x,poly) for poly in polys]
    def jacobian(points):
        columns=[]
        for v,x0,i in zip(points,center,selected,strict=True):
            x=x0*(1+v)
            dif=factors(x)
            g=[arb(0),arb(0)]
            for kind,z,coef in rows:
                d=(-z/x).exp()*dif[kind]/(4*pi)
                for j in range(2):
                    g[j]-=coef[j]*d
            ex=[(-z/x).exp() for _,_,z in field]
            kg=x*sum(a*e for (_,a,_),e in zip(field,ex,strict=True))/(2*pi)
            wg=-x*x*sum(b*e for (b,_,_),e in zip(field,ex,strict=True))/(4*pi)
            columns.append([c[i]*a/(1+v) for a in (*g,kg,wg)]+[x.sqrt()/(2*(1+v))])
        return arb_mat(columns).transpose()


    cols=[values(x) for x in anchors]
    f=[target+sum(w*col[j] for w,col in zip(c,cols,strict=True))
       for j,target in enumerate((*starts,-1,-arb("1e-122")))]
    f.append(sum(x.sqrt() for x in anchors)-sum(source))
    # 对 log x 的单物种角核导数绝对值至多一。第十层引力导数用 0 <= exp(-z/x) <= 1。
    low=anchors[9]*(1+arb(0,width.upper()))
    log_span=-(1-width).log()
    err=[sum(abs(b[j]) for _,_,b in rows)/(4*pi) for j in range(2)]
    err += [low*sum(abs(a) for _,a,_ in field)/(2*pi),low**2*sum(abs(b) for b,_,_ in field)/(4*pi)]
    param=[arb(0,(abs(c[9])*e*log_span).upper()) for e in err]
    param.append(arb(0,(low.sqrt()*log_span/2).upper()))
    f=arb_mat([[value+error] for value,error in zip(f,param,strict=True)])
    pre=jacobian([arb(0)]*5).mid().inv()
    assert not pre.det().contains(0)
    radius=arb("1e-5")
    box=[arb(0,radius.upper())]*5
    defect=arb_mat([[int(i==j) for j in range(5)] for i in range(5)])-pre*jacobian(box)
    offset=max(abs(v).upper() for v in (pre*f).entries())/radius
    slope=max(sum(abs(v).upper() for v in row) for row in defect.tolist())
    print("中心位移占比、导数误差行和",offset.str(13),slope.str(13),flush=True)
    assert offset+slope<1
    output=anchors[:]
    for i in selected:
        output[i]*=1+box[0]
    output[9]=low
    amps=[x.sqrt() for x in output]
    assert all(source[0]<y<source[-1] for y in amps)
    gaps=lib("schur-horn-gravity-holdout")["ky_fan_gaps"](source,amps)
    series=native["weak_series"](56)
    def weak(x):
        return series(x/p0) if x<10**7*p0 else lib_neutral(x,p0,p1)
    lib_neutral=lib("full-neutral-free-cutoff")["weak_column"]
    rest=sum(c[i]*weak(x) for i,x in enumerate(output) if i!=9)
    ds=[rest+c[9]*weak(anchors[9]*(1+s*width)) for s in (-1,1)]
    print("第十层两个参数端点的弱差、最小主序余量",*(v.str(15) for v in ds),min(gaps).str(15),flush=True)
    assert ds[0]>0>ds[1]
    colour_centre=a3+sum(w*col[4] for w,col in zip(c,cols,strict=True))
    colour_lipschitz=sum(abs(b[2]) for _,_,b in rows)/(4*pi)
    variation=sum(abs(c[i]) for i in selected)*(-(1-radius).log())+abs(c[9])*log_span
    colour=colour_centre+arb(0,(colour_lipschitz*variation).upper())
    print("全解盒端点色逆耦合上下界",colour.lower().str(14),colour.upper().str(14),flush=True)
    assert colour>0
    running=lib("full-neutral-running")
    neutral=lib("full-neutral-free-cutoff")
    species,_=neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    err=running["weak_bound"](data,[p0,*output],[arb(0),*c],lib("neutral-axial-top"),lib("neutral-bosonic"))
    for mu in ("1e4","1e9"):
        p=(arb(mu)/arb("1.22089e19"))**2
        centers=[x.mid() for x in output]
        cols=running["direct"](data,centers,c,species,mu,gauge,gauge["angular_rule"](13,37))
        variation=sum(abs(w)*abs((x/mid).log())*min(arb(1),((p-p0)/(4*x.lower())).upper())
                      for w,x,mid in zip(c,output,centers,strict=True))
        for i in range(3):
            r=variation*sum(abs(b[i]) for _,_,b in species)/(4*pi)
            cols[i]+=arb(0,r.upper())
        qq,mix,ac=cols
        det=qq*(mix+arb(0,err.upper()))-mix*mix
        print("能标、QQ、混合、色、行列式",mu,[(v.lower().str(12),v.upper().str(12)) for v in (*cols,det)],flush=True)
        if mu=="1e4":
            assert qq>0 and ac>0 and det>0
            # 在 [8,12] TeV 上逐盒夹住 dR/dlog p。四个角权均为 h=u(1-u) 的一次多项式。
            pl=(arb(8000)/arb("1.22089e19"))**2
            ph=(arb(12000)/arb("1.22089e19"))**2
            angular=((0,6,1),(3,-12,1),(5,12,7),(9,12,11))
            bins=[]
            for m in range(64):
                u,v=arb(m)/128,arb(m+1)/128
                kl,ku=pl*u*(1-u),ph*v*(1-v)
                dh=(v**2-u**2)/2-(v**3-u**3)/3
                bins.append((kl,ku,[2*(a*(v-u)+b*dh)/den for a,b,den in angular]))
            deriv=[arb(0)]*3
            for i,x in enumerate(output):
                for kind,z,b in species:
                    if z>0 and ph/(4*z)<arb("1e-8"):
                        r=arb(0).union(ph/(4*z))
                    else:
                        r=arb(0)
                        for kl,ku,weights in bins:
                            if z.is_zero():
                                frac=arb(1)
                            else:
                                frac=(kl/(z+kl)).union(ku/(z+ku))
                            expo=(-(z+ku)/x.lower()).exp().union((-(z+kl)/x.upper()).exp())
                            r+=weights[kind]*frac*expo
                    for j in range(3):
                        deriv[j]-=c[i]*b[j]*r/(4*pi)
            span=(pl/p).log().union((ph/p).log())
            wide=[v+d*span for v,d in zip(cols,deriv,strict=True)]
            wq,wb,wc=wide
            wd=wq*(wb+arb(0,err.upper()))-wb*wb
            print("8至12 TeV 的QQ、混合、色与行列式",[(v.lower().str(12),v.upper().str(12)) for v in (*wide,wd)],flush=True)
            assert wq>0 and wc>0 and wd>0
        else:
            assert ac<0
    return amps


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
