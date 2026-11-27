"""不设弱差约束，给出自由能标及十九个可酉实现的截止同时匹配两项规范残差和引力的见证。"""
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
    c=[v.real for v in pivots[0]]
    assert sum(c).contains(1)
    moments,heat,(field,)=lib("three-index-free-atoms")["broken_fields"](root,[n+4])
    gauge=lib("all-field-signed-crossing")
    mz=arb("91.1876")
    energy=mz*2**45
    data=gauge["joint_system"](root,energy,number=n+4,order=4,last=8)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,n+4)
    ctx.prec=192
    p0,p1=data["p0"],data["p1"]
    assert p1.overlaps(p0*2**90)
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    species,_=lib("full-neutral-free-cutoff")["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    rows=[(kind,z,(qq-arb(8)*colour/3,mix-colour)) for kind,z,(qq,mix,colour) in species]
    a1,a2,a3=data["initial"]
    starts=(5*a1/3+a2-8*a3/3,a2-a3)
    pi=arb.pi()

    # K=1 单独强制第18个截止较高。其余各层仍允许独立取整个源谱凸盒。
    xmax=ordered[-1]
    tv=sum(abs(w) for w in c)
    beta=sum(a for _,a,_ in field)/(2*pi)
    newton_error=tv*xmax*sum(abs(a)*(1-heat(2,z/xmax)) for _,a,z in field)/(2*pi)
    neg_other=sum(-w for i,w in enumerate(c) if i!=17 and w<0)
    x18=(((1-newton_error)/abs(beta)-neg_other*xmax)/abs(c[17])).lower()
    assert x18>arb(".04907203") and c[17]<0 and beta<0
    pos=sum(w for w in c if w>0)
    pcut=(arb("1e11")/arb("1.22089e19"))**2
    logmax=(pcut/p0).log()
    # 把带电 W 与无质量色规范核配对。正核 H 对 x 单调，且
    # H(xmax)-H(x18)<=(p-p0)/(4*x18)；带质量的差受无截止端点对数矩控制。
    assert [kind for kind,_,_ in rows[:5]]==[0,0,3,1,3]
    for (_,_,b),expected in zip(rows[:5],(0,0,-arb(22)/3,arb(1)/3,arb(88)/3),strict=True):
        assert (b[0]-expected).contains(0)
    mass_w=rows[2][1]
    missing=arb(12)/11*(1+(4*mass_w/p0).log1p())
    total_H=(pos-abs(c[17]))*logmax+abs(c[17])*(pcut-p0)/(4*x18)
    light=22*total_H+pos*(arb(22)/3*missing+logmax/3)
    heavy=tv*sum(abs(b[0])*(pcut-p0)/(4*z) for _,z,b in rows[5:])
    floor=starts[0]-(light+heavy)/(4*pi)
    assert floor>arb("5.01758")
    print("从 M_Z 到 1e11 GeV 的光子残差统一下界",floor.lower().str(12),flush=True)

    selected=(0,15,17,18,16)
    center=list(map(arb,("3.09485665563086e-9","4.03994511350391e-5","0.0815643354325029",
                         "0.368314096285987","0.0818925608720361")))
    fixed=arb("3e-12")
    free_weight=sum(c[i] for i in range(19) if i not in selected)
    rule=gauge["angular_rule"](16,48)
    log=(p1/p0).log()
    def values(x):
        col=[arb(0),arb(0)]
        for kind,z,beta in rows:
            exposure=gauge["angular_exposure"](kind,z,x,p0,p1,log,rule)/(4*pi)
            for j in range(2):
                col[j]-=beta[j]*exposure
        return [*col,*moments(field,heat,x)]
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

    fixed_values=values(fixed)
    cols=[values(x) for x in center]
    f=[target+free_weight*fixed_values[j]+sum(c[i]*col[j] for i,col in zip(selected,cols,strict=True))
       for j,target in enumerate((*starts,-1,-arb("1e-122")))]
    f.append(14*fixed.sqrt()+sum(x.sqrt() for x in center)-sum(source))
    f=arb_mat([[v] for v in f])
    pre=jacobian([arb(0)]*5).mid().inv()
    assert not pre.det().contains(0)
    radius=arb("1e-5")
    box=[arb(0,radius.upper())]*5
    defect=arb_mat([[int(i==j) for j in range(5)] for i in range(5)])-pre*jacobian(box)
    offset=max(abs(v).upper() for v in (pre*f).entries())/radius
    slope=max(sum(abs(v).upper() for v in row) for row in defect.tolist())
    print("中心位移占比、导数误差行和",offset.str(13),slope.str(13),flush=True)
    assert offset+slope<1
    output=[fixed.sqrt()]*19
    for i,x,v in zip(selected,center,box,strict=True):
        output[i]=(x*(1+v)).sqrt()
    assert all(source[0]<y<source[-1] for y in output)
    gaps=lib("schur-horn-gravity-holdout")["ky_fan_gaps"](source,output)
    print("能标、近似中心与最小主序余量",energy.str(17),[v.str(15) for v in center],min(gaps).str(13),flush=True)
    # 见证确定以后读独立弱差；高能端余量对 M>=1e13 GeV 一致。
    lib("neutral-residual-all-cutoffs")["verify"]()
    assert energy>arb("1e13") and all(y*y>p0 for y in output)
    ref=arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    weak_error=tv*arb("2.398e-21")+arb("7.468e-36")*sum(abs(w)/y**2 for w,y in zip(c,output,strict=True))
    weak=ref+arb(0,weak_error.upper())
    assert arb("-.074256577194663365579")<weak<arb("-.074256577194663365511")
    print("同一 M 和输运的独立弱差上下端",weak.lower().str(23),weak.upper().str(23),flush=True)
    return output


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
