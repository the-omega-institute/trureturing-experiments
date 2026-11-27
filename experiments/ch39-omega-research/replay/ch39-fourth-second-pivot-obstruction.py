"""第四壳第二枢轴在原能标下的五条件排除。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root, *, weak_minimum=None):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    n=2_306_888_100_085
    native=lib("native-rh-ray-pvm")
    _,_,pivots,_=next(row for row in native["shell_weights"](root,n+4) if row[0]==6)
    c=[z.real for z in pivots[1]]
    _,ordered=lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=n+1)[2]
    _,heat,(fields,)=lib("three-index-free-atoms")["broken_fields"](root,[n+4])
    gauge=lib("all-field-signed-crossing")
    neutral=lib("full-neutral-free-cutoff")
    data=gauge["joint_system"](root,"1.176397e14",order=4,last=8,number=n+4)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,n+4)
    ctx.prec=192
    p0,p1=data["p0"],data["p1"]
    lo,hi=ordered[0],ordered[-1]
    tv=sum(map(abs,c))
    assert sum(c).contains(1) and all(w<0 or w>0 for w in c)
    assert [i+1 for i,w in enumerate(c) if w>0]==[1,10,16]
    assert min(w for w in c if w>0)>arb("2.51")
    assert max(-w for i,w in enumerate(c) if w<0 and i!=17)<arb("2.913")
    trace=sum(v.sqrt() for v in ordered)
    assert trace<arb("1.186") and abs(c[17])<63
    pi=arb.pi()
    ak=sum(a for _,a,_ in fields)/(2*pi)
    aw=-sum(g for g,_,_ in fields)/(8*pi)
    ek=tv*hi*sum(abs(a)*(1-heat(2,z/hi)) for _,a,z in fields)/(2*pi)
    ew=tv*hi**2*sum(abs(g)*(arb(1)/2-heat(3,z/hi)) for g,_,z in fields)/(4*pi)
    xcut=arb("1e-5")
    d=(1-ek)/abs(ak)-abs(c[17])*xcut
    eps=(arb("1e-122")+ew)/abs(aw)
    assert ak<0<aw and d>arb(".35") and eps<arb(".0003")
    # 第十八层若小，其余负权的一次矩 A 与二次矩 B 满足 Hölder 界。
    # 对正权令 p=Σc x, q²=Σc x²/p。W 控制 B<=p q²+eps，K 给出 A>p+d。
    # q<=Σ正 sqrt(x)，p>=min(c+)q²，负权振幅和至多 T-q。
    # 因而 (p+d)^3 <= (p q²+eps) [2.913(T-q)]²。以下有理盒排尽此必要条件。
    T=arb("1.186")
    nmax=arb("2.913")
    assert nmax*T*T<5
    def exclude(pl,pu,ql,qu):
        if pu < arb("2.51")*ql**2:
            return
        if (pl+d)**3>(pu*qu**2+eps)*(nmax*(T-ql))**2:
            return
        assert max(pu-pl,qu-ql)>arb("1e-4"),(pl,pu,ql,qu)
        if pu-pl>qu-ql:
            m=(pl+pu)/2
            exclude(pl,m,ql,qu)
            exclude(m,pu,ql,qu)
        else:
            m=(ql+qu)/2
            exclude(pl,pu,ql,m)
            exclude(pl,pu,m,qu)
    # 由 A<=2.913T² 与 A>p+d 得 p<5。
    exclude(arb(0),arb(5),arb(0),T)
    print("K,W及振幅迹强制第18层截止平方大于",xcut.str(10),"矩差下限和误差上限",d.str(10),eps.str(10),flush=True)

    series=native["weak_series"](56)
    threshold=8192*p0
    weak_momentum=p1 if weak_minimum is None else (arb(weak_minimum)/arb("1.22089e19"))**2
    assert weak_momentum>=(arb("1e13")/arb("1.22089e19"))**2 or weak_minimum is None
    large=neutral["weak_column"](threshold,p0,weak_momentum)
    whole=large
    assert lo>p0
    for e in range(13):
        step=arb(2)**(e-6)
        for j in range(64):
            ratio=((64+j)*step).union((65+j)*step)
            extra=arb(0,(2*ratio*p0/weak_momentum).upper()) if weak_minimum is not None else arb(0)
            whole=whole.union(series(ratio)+extra)
    if_high=sum(w*(large if i in (0,15) else whole) for i,w in enumerate(c))
    print("两正权层均在阈值以上时的弱差上界",if_high.upper().str(10),flush=True)
    assert if_high<0
    if weak_minimum is not None:
        return c, threshold, xcut, p0, if_high.upper()

    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    species,_=neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    rule=gauge["angular_rule"](10,28)
    log=(p1/p0).log()
    coef=[qq-arb(8)*colour/3 for _,_,(qq,_,colour) in species]
    edges=sorted([lo,threshold,xcut,hi]+[p0*2**k for k in range(4,101,5) if p0*2**k < hi])
    exposed=[[gauge["angular_exposure"](kind,z,x,p0,p1,log,rule)/(4*pi)
              for kind,z,_ in species] for x in edges]
    def hull(columns):
        result=columns[0]
        for v in columns[1:]:
            result=result.union(v)
        return result
    pieces=[-sum(b*l.union(r) for b,l,r in zip(coef,left,right,strict=True))
            for left,right in zip(exposed,exposed[1:])]
    all_x=hull(pieces)
    low=hull(pieces[:edges.index(threshold)])
    high=hull(pieces[edges.index(xcut):])
    a1,a2,a3=data["initial"]
    start=5*a1/3+a2-8*a3/3
    bounds=[start+sum(w*(low if i==small else high if i==17 else all_x)
                      for i,w in enumerate(c)) for small in (0,15)]
    print("第1或第16层小截止时的光子残差下界",*(v.lower().str(12) for v in bounds),flush=True)
    assert all(v>0 for v in bounds)
    return bounds


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
