"""第四壳第三枢轴在原能标下五式相容的区间见证。"""
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
        "0.0059635173156968", "2.11616221952219e-18", "1.39799555210032e-14", "2.0612839852922e-20",
        "4.06913149566367e-31", "7.4286087689067e-05", "2.74842854286246e-18", "0.27194222242332",
        "5.68180702408289e-14", "4.15649725512602e-34", "0.0577119088492197", "0.0869095658763582",
        "1.45752912084902e-23", "5.02124408540481e-18", "1.05205781896682e-28", "6.81068490953953e-10",
        "3.67620555866841e-12", "0.00182320724275718", "2.37954485815224e-21")))
    selected=(17,16,15,7,0)
    center=[anchors[i] for i in selected]
    width=arb('1e-9')
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
    assert arb("29.7103")<colour<arb("30.2330")
    return amps


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
