"""在第三枢轴五式解的同一截止盒上核查中间能标的色、中性动能。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    # 直接复用端点认证返回的振幅盒，固定位与六个调节范围均随之确定。
    sites=[y*y for y in lib("fourth-third-pivot-five")["verify"](root)]
    n=2_306_888_100_085
    native=lib("native-rh-ray-pvm")
    _,_,pivots,_=next(row for row in native["shell_weights"](root,n+4) if row[0]==6)
    weights=[z.real for z in pivots[2]]
    assert sum(weights).contains(1)
    gauge=lib("all-field-signed-crossing")
    energy=arb("1e4")
    data=gauge["joint_system"](root,energy,number=n+4,order=4,last=8)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,n+4)
    ctx.prec=192
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    species,_=lib("full-neutral-free-cutoff")["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    running=lib("full-neutral-running")
    p0,p=data["p0"],data["p1"]
    # 零权辅助位置令通用界对全部十九个位置采用共同无截止端点和有限截止余项。
    err=running["weak_bound"](data,[p0,*sites],[arb(0),*weights],
                              lib("neutral-axial-top"),lib("neutral-bosonic"))
    centers=[x.mid() for x in sites]
    cols=running["direct"](data,centers,weights,species,energy,gauge,gauge["angular_rule"](13,37))
    # 对任一归一化角核，dR/dlog(x) 的绝对值 <= min(1,(p-p0)/(4x))。
    # 这是两个指数之差的上界，使用 1-exp(-h(p-p0)/x)<=h(p-p0)/x。
    movement=sum(abs(w)*abs((x/mid).log())*min(arb(1),((p-p0)/(4*x.lower())).upper())
                 for w,x,mid in zip(weights,sites,centers,strict=True))
    for i in range(3):
        error=movement*sum(abs(beta[i]) for _,_,beta in species)/(4*arb.pi())
        cols[i]+=arb(0,error.upper())
    qq,mix,colour=cols
    pure=mix+arb(0,err.upper())
    det=qq*pure-mix*mix
    print("弱差统一绝对上界",err.upper().str(12),flush=True)
    for name,v in zip(("QQ","QT3","色","T3T3","中性行列式"),(qq,mix,colour,pure,det),strict=True):
        print("10 TeV",name,v.lower().str(13),v.upper().str(13),flush=True)
    assert qq>0 and pure<0 and colour<0 and det<0
    a1,a2,a3=data["initial"]
    assert min(a1,a2,a3)>0
    return colour,pure,det


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
