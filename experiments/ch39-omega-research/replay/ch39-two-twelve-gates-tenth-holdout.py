"""按预登记协议比较同门序双九式见证传至第十壳的三个端点。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    first = lib("twelve-independent-gates-nine")
    second = lib("second-twelve-gates-nine")
    transport = lib("two-shell-ordered-bond-gates")["transport"]
    layers = tuple((j,(edge,)) for j,edge in enumerate(first["PAIRS"]))
    boxes = [item["experiment"](root,fifth,columns,levels)[0] for item in (first,second)]
    _,_,_,basis,_,weights,heat,field = lib("fifth-seventh-shortest-site-gates")["geometry"](
        root,fifth,columns,levels,number=2_306_888_100_095)
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    weak = lib("full-neutral-free-cutoff")["weak_column"]
    p0 = (arb("91.1876")/arb("1.22089e19"))**2
    p1 = (arb("1.176397e14")/arb("1.22089e19"))**2
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    left = basis.conjugate().transpose()
    source = fifth*columns
    out = []
    assert sum(weights).contains(1)
    for k,(model,box) in enumerate(zip((first,second),boxes,strict=True)):
        def projected(point):
            return [left*g for g in transport(source,point,parameters=model["PARAMETERS"],
                                              active=model["ACTIVE"],layers=layers)]
        center = [x.mid() for x in box]
        radii = [(x-c).abs_upper() for x,c in zip(box,center,strict=True)]
        point,*center_jac = projected(center)
        block,*derivs = projected(box)
        f,_ = gravity(point,center_jac,levels,weights,heat,field,True)
        _,j = gravity(block,derivs,levels,weights,heat,field,False)
        grav = [f[i,0]+target+sum(j[i,a]*(x-c) for a,(x,c) in enumerate(zip(box,center,strict=True)))
                for i,target in enumerate((arb(1),arb("1e-122")))]
        def lower(i,c):
            movement = sum(r*abs(d[i,c]).upper() for r,d in zip(radii,derivs,strict=True))
            return max(arb(0),(abs(point[i,c]).lower()-movement).lower())
        low,heavy,top = levels
        amps = [low+(heavy-low)*lower(i,0)**2+(top-low)*lower(i,1)**2 for i in range(19)]
        assert min(a.lower() for a in amps)>p0.sqrt()
        discrepancy = ref+sum(w*(weak(a.lower()**2,p0,p1)-ref) for a,w in zip(amps,weights,strict=True))
        out.append([*grav,discrepancy])
        print("见证、牛顿、真空、弱端点",k+1,
              [(v.lower().str(20),v.upper().str(20)) for v in out[-1]],flush=True)
    delta = [b-a for a,b in zip(*out,strict=True)]
    print("第二组减第一组",[(d.lower().str(20),d.upper().str(20)) for d in delta],flush=True)
    bands = ((("1.01842","1.02377"),("-.002557","-.002123")),
             (("1.03959","1.04895"),("-.000581","-.0000926")))
    for values,bounds in zip(out,bands,strict=True):
        for value,(lo,hi) in zip(values[:2],bounds,strict=True):
            assert arb(lo)<value<arb(hi)
        assert arb("-.074256577194663371")<values[2]<arb("-.074256577194663360")
    assert arb(".01582538")<delta[0]<arb(".03052422")
    assert arb(".00154322")<delta[1]<arb(".00246405")
    assert delta[2].contains(0) and abs(delta[2])<arb("9.337e-18")
    return out,delta


def verify(root):
    original = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    return experiment(root,*original["verify"](root,return_operator="spectral"))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
