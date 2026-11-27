"""第四壳三低调整盒的端点色逆耦合与连续路径上的动能符号。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    n=2_306_888_100_085
    native=lib("native-rh-ray-pvm")
    _,_,pivots,_=next(row for row in native["shell_weights"](root,n+4) if row[0]==6)
    c=[z.real for z in pivots[2]]
    gauge=lib("all-field-signed-crossing")
    energy=arb("1.176397e14")
    data=gauge["joint_system"](root,energy,number=n+4,order=4,last=8)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,n+4)
    ctx.prec=192
    box=lib("fourth-third-pivot-three-adjusters")["local_box"]()
    mid=[arb(0)]*19
    radius=[box["high_log_radius"]]*19
    for i,x,r in zip(box["layers"],box["centers"],box["radii"],strict=True):
        mid[i],radius[i]=x,r
    high=iter(box["high_centers"])
    for i in range(19):
        if i not in box["layers"]:
            mid[i]=next(high)
    assert next(high,None) is None
    sites=[x*arb(0,r.upper()).exp() for x,r in zip(mid,radius,strict=True)]
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    species,_=lib("full-neutral-free-cutoff")["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    running=lib("full-neutral-running")
    p0=data["p0"]
    err=running["weak_bound"](data,[p0,*sites],[arb(0),*c],lib("neutral-axial-top"),lib("neutral-bosonic"))
    rule=gauge["angular_rule"](17,47)
    initial1,initial2,initial3=data["initial"]
    assert min(initial1,initial2,initial3)>0
    print("减法点的色逆耦合及中性行列式",initial3.str(14),(5*initial1*initial2/3).str(14),flush=True)
    for mu in (energy,arb("1e4")):
        p=(mu/arb("1.22089e19"))**2
        cols=running["direct"](data,mid,c,species,mu,gauge,rule)
        # 每个物种归一化正角权给出 |dR/dlog x|<=min(1,(p-p0)/(4x))，
        # 在高端十六个截止的逐坐标完整对数盒上累加误差。
        movement=sum(abs(w)*r*min(arb(1),((p-p0)/(4*x.lower())).upper())
                     for w,r,x in zip(c,radius,sites,strict=True))
        for j in range(3):
            delta=movement*sum(abs(beta[j]) for _,_,beta in species)/(4*arb.pi())
            cols[j]+=arb(0,delta.upper())
        qq,mix,colour=cols
        pure=mix+arb(0,err.upper())
        det=qq*pure-mix*mix
        print("能标、QQ、混合、色、纯弱、行列式",mu.str(9),[(v.lower().str(14),v.upper().str(14))
              for v in (*cols,pure,det)],flush=True)
        if mu>arb("1e13"):
            assert qq>0 and colour>0 and det>0
        else:
            assert qq>0 and colour<0 and det<0
    return True


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
