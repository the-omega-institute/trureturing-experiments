"""同一完整夸克第四步精确复现全程正性谱的六个强权截止。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root, *, weak_zero=False):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    previous = lib("fifth-zero-strong-flavor")
    data,species,_,_,boxes,_ = lib("fifth-third-pivot-three-high-positive")["endpoints"](root)
    rows = (0,9,15,16,17,18)
    def block_target(root, lib, blocks):
        return lib("fifth-free-frame-blocks")["solve"](
            root,lib,blocks,previous["NAME"],previous["COEFFICIENTS"],previous["FRAME"])

    if weak_zero:
        from inspect import getclosurevars
        number = 2_306_888_100_090
        weights=next(z[2][2] for z in lib("native-rh-ray-pvm")["shell_weights"](root,number) if z[0]==6)
        q0,q1=data["p0"],data["p1"]
        upper=arb("1e-5")
        weak=lib("full-neutral-free-cutoff")["weak_column"]
        rest=arb(0)
        for i,w in enumerate(weights):
            if i in rows:
                continue
            rest += w.real*(weak(upper,q0,q1)-weak(boxes[i],q0,q1))
        series=lib("native-rh-ray-pvm")["weak_series"](56)
        coefficients=getclosurevars(series).nonlocals["coefficients"]
        width=arb("1e-9")
        for sign in (-1,1):
            shift=sign*width
            inverse=q0/boxes[0]
            correction=sum(a*inverse**n*(-n*shift).expm1()
                           for n,a in enumerate(coefficients,1))
            tail=arb(0,(arb("1e-19")+4*boxes[0]*width.exp()/q1).upper())
            drift=rest+weights[0].real*(correction+tail)
            print("第一层对数变化、弱差端点严格包围",sign,drift.str(10),flush=True)
            assert drift*sign<0
        boxes=list(boxes)
        boxes[0]*=arb(0,width.upper()).exp()
    sites = lib("fifth-flavor-cutoff-chart")["solve"](
        root, "ch39-fifth-six-positive-flavor", (0,0,0,1),previous["DIRECTIONS"],frame=block_target,
        rows=rows,target_sites=[boxes[i] for i in rows])

    if weak_zero:
        assert min(x.lower() for i,x in enumerate(sites) if i not in rows)>upper
        return sites
    # 原见证的弱差精确为零，两者强权截止相同，只剩十三个弱权层的差。
    number = 2_306_888_100_090
    weights=next(z[2][2] for z in lib("native-rh-ray-pvm")["shell_weights"](root,number) if z[0]==6)
    p0=(arb("91.1876")/arb("1.22089e19"))**2
    p1=(arb("1.176397e14")/arb("1.22089e19"))**2
    weak=lib("full-neutral-free-cutoff")["weak_column"]
    weights=[w.real for w in weights]
    difference=sum(weights[i]*(weak(sites[i],p0,p1)-weak(boxes[i],p0,p1))
                   for i in range(19) if i not in rows)
    print("同六个强权截止时的端点弱差区间", difference.str(12),flush=True)
    assert abs(difference) < arb("3e-12")
    for i in rows:
        sites[i] = boxes[i]
    moments,heat,(field,) = lib("three-index-free-atoms")["broken_fields"](root,[number])
    thermal = [moments(field,heat,x) for x in sites]
    values = [sum(c*r[j] for c,r in zip(weights,thermal,strict=True)) for j in range(2)]
    gauge=lib("all-field-signed-crossing")
    rule=gauge["angular_rule"](17,46)
    log_ratio=(p1/p0).log()
    end=list(data["initial"])
    end=[5*end[0]/3+end[1],end[1],end[2]]
    for i,(x,c) in enumerate(zip(sites,weights,strict=True)):
        point=x.mid()
        delta=abs((x/point).log()).upper()
        for kind,z,b in species:
            exposure=gauge["angular_exposure"](kind,z,point,p0,p1,log_ratio,rule)+arb(0,delta)
            for j in range(3):
                end[j]-=c*b[j]*exposure/(4*arb.pi())
    rp,rm=end[0]-8*end[2]/3,end[1]-end[2]
    print("六层匹配后的牛顿项、真空项、端点色读数与两个规范残差",
          *[x.str(10) for x in [*values,end[2],rp,rm]],flush=True)
    assert arb(".7")<values[0]<arb(".85") and arb("-.45")<values[1]<arb("-.30")
    assert end[2]<-60 and rp>100 and rm>50
    return sites


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
