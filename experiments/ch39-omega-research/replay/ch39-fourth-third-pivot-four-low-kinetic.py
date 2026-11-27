"""四枚低截止在第三枢轴上同时给出三个端点等式及 5.6 TeV 附近的正动能。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    n=2_306_888_100_085
    native=lib("native-rh-ray-pvm")
    _,_,pivots,_=next(row for row in native["shell_weights"](root,n+4) if row[0]==6)
    c=[z.real for z in pivots[2]]
    _,source_x=lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=n+1)[2]
    source=[x.sqrt() for x in source_x]
    energy=arb("1.176397e14")
    gauge=lib("all-field-signed-crossing")
    data=gauge["joint_system"](root,energy,number=n+4,order=4,last=8)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,n+4)
    ctx.prec=192
    neutral=lib("full-neutral-free-cutoff")
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    species,_=neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    pi=arb.pi()
    rows=[(k,z,[(qq-8*col/arb(3))/(4*pi),(mix-col)/(4*pi)]) for k,z,(qq,mix,col) in species]
    p0,p1=data["p0"],data["p1"]
    a1,a2,a3=data["initial"]
    # 四个低层，先固定第 17 层，再用第 10 层弱差参数和第 1、18 层规范变量。
    layers=(0,9,16,17)
    anchor=list(map(arb,("1.946148411827196e-15","3.56255408181904e-34",
                         "3.308558011555144e-33","4.222671821839275e-10")))
    radii=list(map(arb,("1e-5","5e-9","0","1e-5")))
    lows=[x*arb(0,r.upper()).exp() for x,r in zip(anchor,radii,strict=True)]
    def fill(xs):
        return ((sum(source)-sum(x.sqrt() for x in xs))/15)**2
    high_mid,high=fill(anchor),fill(lows)
    wh=sum(c[i] for i in range(19) if i not in layers)
    assert all(x<arb("1e-6") for x in lows) and high>arb("1e-6")
    rule=gauge["angular_rule"](19,51)
    def gauge_pair(x):
        col=[arb(0),arb(0)]
        for k,z,b in rows:
            h=gauge["angular_exposure"](k,z,x,p0,p1,(p1/p0).log(),rule)
            for j in range(2):
                col[j]-=b[j]*h
        return col
    fm=[5*a1/3+a2-8*a3/3,a2-a3]
    for x,w in zip([*anchor,high_mid],[*[c[i] for i in layers],wh],strict=True):
        col=gauge_pair(x)
        for j in range(2):
            fm[j]+=w*col[j]
    order=24
    beta=[arb(factorial(j))**2/factorial(2*j+1) for j in range(order+2)]
    angular=((0,6,1),(3,-12,1),(5,12,7),(9,12,11))
    polys=[[(a*beta[j]+b*beta[j+1])*(-1)**j/(d*factorial(j)) for j in range(order+1)]
           for a,b,d in angular]
    def average(t,poly):
        if t>100:
            return arb(0).union(12/t)
        assert t<1
        v=arb(0)
        for coeff in reversed(poly):
            v=v*t+coeff
        return v+arb(0,((t/4)**(order+1)/factorial(order+1)).upper())
    def derivative(x):
        def avg(t,poly):
            return average(t.lower(),poly).union(average(t.upper(),poly))
        factors=[avg(p0/x,poly)-avg(p1/x,poly) for poly in polys]
        return [-sum(b[j]*(-z/x).exp()*factors[k] for k,z,b in rows) for j in range(2)]
    def jacobian(xs):
        return arb_mat([[c[layers[i]]*v for v in derivative(xs[i])] for i in (0,3)]).transpose()
    pre=jacobian(anchor).mid().inv().mid()
    assert not pre.det().contains(0)
    disp=pre*arb_mat([[v] for v in fm])
    defect=arb_mat([[1,0],[0,1]])-pre*jacobian(lows)
    tgrad=pre*arb_mat([[c[9]*v] for v in derivative(lows[1])])
    dh=abs((high/high_mid).log()).upper()
    for j in range(2):
        eff=sum(abs(sum(pre[j,k]*b[k] for k in range(2))) for _,_,b in rows)
        maxmove=abs(disp[j,0])+abs(tgrad[j,0])*radii[1]+abs(wh)*eff*(p1-p0)*dh/(4*high.lower())
        ratio=maxmove/radii[0]+sum(abs(defect[j,k]) for k in range(2))
        print("规范更新的自映行占比",j+1,ratio.upper().str(12),flush=True)
        assert ratio<1
    weak=native["weak_series"](56)
    rest=wh*neutral["weak_column"](high,p0,p1)+c[16]*weak(lows[2]/p0)
    rest+=sum(c[layers[i]]*neutral["weak_column"](lows[i],p0,p1) for i in (0,3))
    ds=[rest+c[9]*weak(anchor[1]*(s*radii[1]).exp()/p0) for s in (-1,1)]
    assert ds[0]>0>ds[1]
    positions=[high]*19
    for i,x in zip(layers,lows,strict=True):
        positions[i]=x
    gaps=lib("schur-horn-gravity-holdout")["ky_fan_gaps"](source,[x.sqrt() for x in positions])
    print("第十层弱差两端、公共高截止与最小主序余量",[(v.lower().str(13),v.upper().str(13))
          for v in [*ds,high,min(gaps)]],flush=True)
    running=lib("full-neutral-running")
    # 十五个相等高层先合并，弱差无截止项用归一化，有限截止余量按绝对值界。
    weights=[*[c[i] for i in layers],wh]
    sites=[*lows,high]
    err=running["weak_bound"](data,[p0,*sites],[arb(0),*weights],lib("neutral-axial-top"),lib("neutral-bosonic"))
    for mu in (arb(5600),arb(1000),energy):
        m=[x.mid() for x in sites]
        cols=running["direct"](data,m,weights,species,mu,gauge,rule)
        p=(mu/arb("1.22089e19"))**2
        variation=sum(abs(w)*abs((x/center).log())*min(arb(1),((p-p0)/(4*x.lower())).upper())
                      for w,x,center in zip(weights,sites,m,strict=True))
        for j in range(3):
            dp=variation*sum(abs(b[j]) for _,_,b in species)/(4*pi)
            cols[j]+=arb(0,dp.upper())
        if mu==5600:
            # 直接控制外探动量变化。轻重物种各对 dR/dlog p 给出 min(1,p/(4z))。
            pu=(arb(5610)/arb("1.22089e19"))**2
            span=-2*(arb(5590)/5600).log()
            tv=sum(abs(w) for w in weights)
            for j in range(3):
                e=span*tv*sum(abs(b[j])*(1 if z.is_zero() else min(arb(1),(pu/(4*z)).upper()))
                              for _,z,b in species)/(4*pi)
                cols[j]+=arb(0,e.upper())
        qq,mix,colour=cols
        det=qq*(mix+arb(0,err.upper()))-mix*mix
        print("探测能标、QQ、混合、色、行列式",mu.str(10),
              [(v.lower().str(13),v.upper().str(13)) for v in (*cols,det)],flush=True)
        if mu==1000:
            assert colour<0
        else:
            assert colour>0 and qq>0 and det>0


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
