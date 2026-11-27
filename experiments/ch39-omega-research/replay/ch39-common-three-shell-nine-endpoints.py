"""同一个保留源谱的荷六算子在第七、八、九壳实现六个引力值与三个弱端点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


ACTIVE = (53,17,61,51,32,62)
PARAMETERS = (
    4.5, 3.9375, -.1875, -1.3125, -.8125, -3, -.0625, -1.0625,
    .1875, .5, 3.5625, 4.6875, -1.375, -2.375, 3.5, 2, 4.1875, "1.1586908516988268",
    -1.625, 4.625, 3.375, 1.875, -2.8125, 1.3125, -4.1875, -2.125, .5, -2.625, -.1875,
    -2.1875, -1.9375, 6.375, "-4.4980578011411625", -2.875, 3.0625, -4.6875, -2.5, 1.375,
    4.6875, -1.125, -3.1875, -1.5, -1.875, 1.5, 0, 3.875, 2.5, 1.375, -2.25, .25,
    -2.8125, "5.406014641710984", -.875, "-3.9086480728137953", -2.125, -.1875, 1.625,
    -2.9375, -.1875, -2.5, -1.1875, "-.12677947022899566", "-2.0880943420322695", -4.5,
    1.1875, 5.3125, -1.5625, 2, 5.25, 6, 1.1875, .3125, -2.3125, -5.1875, -1.4375, -2.625)
Q_CENTERS = ("954.616970399", "895.563345345", "835.615125179")


def adj(a):
    return a.conjugate().transpose()


def inner(a,b):
    return (adj(a)*b)[0,0]


def normalize(col, tangent):
    norm = inner(col,col).real.sqrt()
    assert norm > 0
    unit = col*(1/norm)
    return unit, [(d-unit*inner(unit,d).real)*(1/norm) for d in tangent]


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    _, squared = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root, refine_node=2_306_888_100_086)[2]
    levels = [squared[i].sqrt() for i in (0,17,18)]
    low, heavy, _ = levels
    native, graph = lib("native-rh-ray-pvm"), lib("three-ray-chiral-boundary")
    shells = []
    for n in range(2_306_888_100_092,2_306_888_100_095):
        nodes, basis = lib("adjacent-pvm-transport")["spectrum"](graph,6,n)
        _, other, pivots, _ = next(r for r in native["shell_weights"](root,n) if r[0]==6)
        assert all(a.overlaps(b) for a,b in zip(nodes,other,strict=True))
        _, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root,[n])
        shells.append((adj(basis), [w.real for w in pivots[2]], heat, field))
    ctx.prec = 256
    identity19 = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    rows = acb_mat([[left[15,i] for i in range(19)] for left,_,_,_ in shells])
    gram = rows*adj(rows)
    assert not gram.det().contains(0)
    dual = adj(rows)*gram.inv()
    proj = identity19-dual*rows
    source_parts = [(k//38,(k%38)//2,k%2) for k in ACTIVE]
    ds = [[acb_mat(19,1) for _ in ACTIVE] for _ in range(2)]
    for j,(imag,i,h) in enumerate(source_parts):
        ds[h][j] = acb_mat([[proj[t,i]*(acb(0,1) if imag else 1)] for t in range(19)])
    planck, mz = arb("1.22089e19"), arb("91.1876")
    gravity = lib("shortest-common-complex-gate")["gravity_system"]

    def frame(point, q):
        raw = list(map(arb,PARAMETERS))
        for k,value in zip(ACTIVE,point,strict=True):
            raw[k] = value
        a,b = [proj*acb_mat([[acb(raw[2*i+h], raw[38+2*i+h])] for i in range(19)]) for h in range(2)]
        u,du = normalize(a,ds[0])
        overlap = inner(u,b)
        db = [d-e*overlap-u*(inner(e,b)+inner(u,d)) for e,d in zip(du,ds[1],strict=True)]
        v,dv = normalize(b-u*overlap,db)
        # 固定三条指定谱线上第一重列的实投影，第二重列在共同正交补内。
        amounts = acb_mat([[((energy/planck-low)/(heavy-low)).sqrt()] for energy in q])
        lift = dual*amounts
        scale = (1-inner(lift,lift).real).sqrt()
        assert scale > 0
        def pair(a,b):
            return acb_mat([[a[i,0],b[i,0]] for i in range(19)])
        return [pair(lift+u*scale,v)] + [pair(e*scale,d) for e,d in zip(du,dv,strict=True)]

    def system(point, q, values):
        mats = frame(point,q)
        fs, js = [], []
        for left,weights,heat,field in shells:
            p,*derivatives = [left*m for m in mats]
            f,j = gravity(p,derivatives,levels,weights,heat,field,values)
            fs.extend(f.tolist())
            js.extend(j.tolist())
        return arb_mat(fs), arb_mat(js)

    center = [arb(PARAMETERS[k]) for k in ACTIVE]
    q0 = list(map(arb,Q_CENTERS))
    energy_width = arb(".000001")
    qs = [q+arb(0,energy_width) for q in q0]
    for _ in range(3):
        f,j = system(center,q0,True)
        step = j.mid().solve(f.mid()).entries()
        center = [(x-d).mid() for x,d in zip(center,step,strict=True)]
    f,j = system(center,qs,True)
    pre = j.mid().inv().mid()
    assert not pre.det().contains(0)
    radius = arb("1e-7")
    box = [c+arb(0,radius) for c in center]
    _,jbox = system(box,qs,False)
    identity = arb_mat([[int(i==j) for j in range(6)] for i in range(6)])
    ratios = [v.abs_upper()/radius+sum(x.abs_upper() for x in row)
              for v,row in zip((pre*f).entries(),(identity-pre*jbox).tolist(),strict=True)]
    print("三弱能标盒内六引力式的统一自映比", [s.str(10) for s in ratios],flush=True)
    assert max(ratios) < arb(".000253")
    source = frame(box,qs)[0]
    weak = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    ref = arb("-.074256577194663365545")
    p0,p1 = (mz/planck)**2,(arb("1.176397e14")/planck)**2
    for q,(left,weights,_,_) in zip(q0,shells,strict=True):
        p = left*source
        cutoff = [low+sum((levels[j+1]-low)*abs(p[i,j])**2 for j in range(2)) for i in range(19)]
        minimum = min(x.lower() for i,x in enumerate(cutoff) if i != 15)
        assert minimum > arb(".0003")
        assert sum(weights).contains(1)
        tail = sum(weights[i]*(large(minimum**2,p0,p1)-ref) for i in range(19) if i != 15)
        faces = [ref+weights[15]*(weak(((q+sign*energy_width)/mz)**2)-ref)+tail for sign in (-1,1)]
        print("第十六线能标中心、其余振幅下界及该壳两个面的弱差",
              q.str(16),minimum.str(12),faces[0].str(10),faces[1].str(10),flush=True)
        assert faces[0] > 0 > faces[1]
    return source, levels


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
