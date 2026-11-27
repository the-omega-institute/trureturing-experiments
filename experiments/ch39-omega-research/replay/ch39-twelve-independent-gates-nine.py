"""平方弦长至多 52 的十二条边各取一次独立二站酉门，三壳的九个端点。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


PAIRS = ((6,11),(7,12),(8,18),(0,2),(1,3),(4,5),
         (6,9),(7,15),(8,16),(6,13),(7,10),(8,17))
LINES = (0,15,0)
ENERGIES = ("522.398135900","895.563345345","544.705683199")
PARAMETERS = (
    "1.25", "-2.3125", "-2.8726920441021", "2.6875",
    "-1.7103926750705", "-2.875", "1.6995056579175", "-0.5",
    "0.625", "2.3125", "-0.5", "1.9375",
    "-2.680818063816", "0.9570626810273", "1.5625", "-2.8307560621693",
    "-0.49352311143283", "-0.66221627224785", "-0.875", "-3",
    "2.1608093428322", "-2.0932731671486", "-2.125", "-0",
    "-1.8811840048849", "-2", "-0.16300248959085", "-1.1875",
    "0.3904006855994", "2", "0.8125", "-1.75",
    "-1.3809017156157", "2.9375", "-2.125", "-2.5",
    "0.59207229125262", "0.4375", "-1.6875", "2.9375",
    "2.6759536211169", "0.0625", "-0.125", "2.875",
    "-1.8703404505959", "-2.1875", "-1.3125", "0.69765039319415",
)
ACTIVE = (2, 4, 6, 12, 13, 15, 16, 17, 20, 21, 24, 26, 28, 32, 36, 40, 44, 47)


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    ordered = lib("two-shell-ordered-bond-gates")
    coords = lib("three-ray-chiral-boundary")["geo"]["sites"]
    short_edges = {(a,b) for a in range(19) for b in range(a+1,19)
                   if sum((u-v)**2 for u,v in zip(coords[a],coords[b],strict=True))<=52}
    assert len(PAIRS)==len(short_edges)==12 and {tuple(sorted(e)) for e in PAIRS}==short_edges
    layers = tuple((j,(e,)) for j,e in enumerate(PAIRS))
    source = fifth*columns
    geometry = lib("fifth-seventh-shortest-site-gates")["geometry"]
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    shells = []
    for n in range(2_306_888_100_092,2_306_888_100_095):
        _,_,_,basis,_,weights,heat,field = geometry(root,fifth,columns,levels,number=n)
        shells.append((basis.conjugate().transpose(),weights,heat,field))
    mp,mz = arb("1.22089e19"),arb("91.1876")
    low,heavy,top = levels
    q = [arb(v) for v in ENERGIES]
    width = arb(".000001")
    energies = [x+arb(0,width) for x in q]
    def evaluate(point,qs,values):
        mats = ordered["transport"](source,point,parameters=PARAMETERS,active=ACTIVE,layers=layers)
        fs,js,fg,jg = [],[],[],[]
        for (basis,weights,heat,field),line,energy in zip(shells,LINES,qs,strict=True):
            proj,*ds = [basis*m for m in mats]
            b = ((energy/mp-low)/(heavy-low)).sqrt()
            for c in (0,1):
                residual = proj[line,c]-(b if c==0 else 0)
                fs.extend((residual.real,residual.imag))
                js.extend([[v[line,c].real for v in ds],[v[line,c].imag for v in ds]])
            f,j = gravity(proj,ds,levels,weights,heat,field,values)
            fg.extend(f.entries())
            jg.extend(j.tolist())
        return arb_mat([[v] for v in fs+fg]),arb_mat(js+jg)
    center = [arb(PARAMETERS[j]) for j in ACTIVE]
    for _ in range(3):
        f,j = evaluate(center,q,True)
        step = j.mid().solve(f.mid()).entries()
        center = [(x-d).mid() for x,d in zip(center,step,strict=True)]
    f,j = evaluate(center,energies,True)
    pre = j.mid().inv().mid()
    assert not pre.det().contains(0)
    offsets = [v.abs_upper() for v in (pre*f).entries()]
    radius = arb(".00045")
    box = [c+arb(0,radius) for c in center]
    curv = lib("three-shell-gravity-radius56")
    value,grad,second,keys = curv["hessian_columns"](source,layers,box,ordered["block"],
                                                     parameters=PARAMETERS,active=ACTIVE)
    grav_h,_ = curv["gravity_hessian"](value,grad,second,keys,shells,levels)
    projection_h = []
    for (basis,_,_,_),i in zip(shells,LINES,strict=True):
        projected = [basis*m for m in second]
        for c in range(2):
            projection_h += [[p[i,c].real for p in projected],[p[i,c].imag for p in projected]]
    h = pre*arb_mat(projection_h+grav_h.tolist())
    eye = arb_mat([[int(i==k) for k in ACTIVE] for i in ACTIVE])
    center_defect = [sum(z.abs_upper() for z in row) for row in (eye-pre*j).tolist()]
    second_norm = [sum(z.abs_upper()*(1 if a==b else 2) for z,(a,b) in zip(row,keys,strict=True))
                   for row in h.tolist()]
    ratio = max(o/radius+e+radius*g/2 for o,e,g in zip(offsets,center_defect,second_norm,strict=True))
    contraction = max(e+radius*g for e,g in zip(center_defect,second_norm,strict=True))
    print("中心偏移、十八维自映比与压缩率",
          *[v.str(12) for v in (max(offsets),ratio,contraction)],flush=True)
    assert ratio < arb(".764") and contraction < arb(".325")
    center_mats = ordered["transport"](source,center,parameters=PARAMETERS,active=ACTIVE,layers=layers)
    native = lib("native-rh-ray-pvm")
    finite = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    for (basis,weights,_,_),i,mid in zip(shells,LINES,q,strict=True):
        centered = basis*center_mats[0]
        directions = [basis*g for g in grad]
        def abs_lower(k,c):
            motion = radius*sum(abs(d[k,c]).upper() for d in directions)
            return max(arb(0),(abs(centered[k,c]).lower()-motion).lower())
        smallest = min((low+(heavy-low)*abs_lower(k,0)**2+(top-low)*abs_lower(k,1)**2).lower()
                       for k in range(19) if k!=i)
        assert smallest > p0.sqrt()
        tail = (large(smallest**2,p0,p1)-ref).abs_upper()*sum(abs(w).upper() for k,w in enumerate(weights) if k!=i)
        faces = [ref+weights[i]*(finite(((mid+sign*width)/mz)**2)-ref)+arb(0,tail) for sign in (-1,1)]
        print("低线、其余振幅下界、能标盒的两个弱差面",i+1,*[v.str(12) for v in (smallest,*faces)],flush=True)
        assert faces[0]>0>faces[1]
    return box,energies


def verify(root):
    original = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    return experiment(root,*original["verify"](root,return_operator="spectral"))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
