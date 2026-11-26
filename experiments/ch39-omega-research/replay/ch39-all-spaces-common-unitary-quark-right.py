"""在右侧实际壳覆盖全部四站几何与共同弱基，排除同一十九维酉运输。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


CAPS = (".947555", ".293784")
TRACE = ".0119"
CHARACTER = ".00075"
ORBITS = ((0, 3, 4), (1, 2, 5), (6, 7, 8), (9, 15, 16),
          (10, 13, 17), (11, 12, 18), (14,))


def certify_observables(root):
    common = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    earlier = run_path(str(root / "ch39-common-unitary-quark-right-trace.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    spectrum = run_path(str(root / "ch39-adjacent-pvm-transport.py"))["spectrum"]
    adj, weighted, eigs = common["adj"], earlier["weighted"], earlier["eigs"]
    source = 2_306_888_100_086
    _, label, _ = common["verify_nodes"](root, ((source,".0294"),), ".1", (),
                                        ((6,7,10,13), (6,7,10,12)))
    ctx.prec = 384
    spaces = [acb_mat([[1/arb(len(g)).sqrt() if i in g else 0 for g in ORBITS]
                       for i in range(19)])]
    sites = shell["geo"]["sites"]
    for sign in (-1, 1):
        omega = acb(-arb(1)/2, sign*arb(3).sqrt()/2)
        v = acb_mat(19, 6)
        for j, orbit in enumerate(ORBITS[:-1]):
            i = orbit[0]
            for n in range(3):
                v[i,j] = omega**n/arb(3).sqrt()
                x,y,z = sites[i]
                i = sites.index((y,z,x))
            assert i == orbit[0]
        spaces.append(v)
    assert ORBITS[2] == (6,7,8) and ORBITS[6] == (14,)
    operators, obs = [], []
    for charge, cap in zip((-4,2), CAPS, strict=True):
        _, before = spectrum(shell, charge, source)
        _, after = spectrum(shell, charge, source+1)
        w = label[source,charge]
        assert all(0 < z < arb(cap) for z in w)
        operators.append(weighted(before,w))
        measured = []
        for k, basis in enumerate(spaces):
            wave = adj(basis)*after
            # 极点和对称三射线给平凡块迹；另两个坐标为三射线的一维 Fourier 投影。
            indices = (2,6) if k == 0 else (2,)
            readout = [sum(abs(wave[l,j])**2 for l in indices) for j in range(19)]
            measured.append(weighted(wave, readout))
        obs.append(measured)
    h = operators[0]+15*operators[1]
    lower = []
    for k, basis in enumerate(spaces):
        block = adj(basis)*h*basis
        if k == 0:
            lower.extend(z.lower() for z in eigs(block))
        else:
            edge = min((block[i,i].real-sum(abs(block[i,j]) for j in range(6) if j != i)).lower()
                       for i in range(6))
            lower.extend([edge]*6)
    lower.sort()
    assert lower[0] > 0
    gaps = []
    for k, (a,b) in enumerate(zip(*obs, strict=True)):
        earlier["positive"](a)
        earlier["positive"](b)
        g = (a.inv()+arb(15)/25*b.inv()).inv()
        vals = eigs(g)
        assert vals[0] > 0
        bound = sum(x*y.lower() for x,y in zip(lower[:len(vals)], reversed(vals), strict=True))
        threshold = arb(CHARACTER if k else TRACE)
        assert bound > threshold
        gaps.append((bound-threshold)/(arb(CAPS[0])*(1 if k else 2)))
        print("表示块、秩、两荷加权下界",k,len(vals),bound.str(16),flush=True)
    gap = min(v.lower() for v in gaps)
    assert gap > arb(".000319")
    print("全部空间上双荷输运的最小射影算符范数距离下界",gap.str(15),flush=True)


def cover_orientations(model, split, poly, blocks, h, chart, rows):
    stack = [(poly,0)]
    while stack:
        p, depth = stack.pop()
        inv = model["box_invariants"](p, blocks, h, chart)
        if inv is False or (inv and all(model["ray_excluded"](row,inv) for row in rows)):
            continue
        if depth == 8:
            return False
        pieces = split(model,p,h)
        if not pieces:
            return False
        stack.extend((q,depth+1) for q in pieces if q)
    return True


def exclude(model, split, blocks, chart, t, rows):
    result = model["polygon"](blocks,chart,t)
    if result is None:
        return False
    poly, h = result
    if not poly:
        return True
    remote = 1-chart
    _,s,z,_,w,_ = blocks[remote]

    def cut(a,b,c):
        nonlocal poly
        poly = model["cut"](poly,*model["coeff"](a,b,c))

    min_u, min_v = [], []
    for j,(d0,s0,_,k0,w0,_) in enumerate(blocks):
        cap = arb(CAPS[j])
        if s0 > cap or d0 > cap:
            return True
        if not (s0 > 0 and s0 < cap and d0 < cap):
            return False
        min_u.extend((w0/s0, -w0/(cap-s0)))
        min_v.extend((-w0/s0, w0/(cap-s0)))
        coeff = abs(k0)**2/(cap-d0)
        if j == remote:
            cut(coeff,arb(1),cap-coeff)
        else:
            cut(coeff+h[1],h[2],cap-coeff-h[0])
    mean = blocks[0][1]+25*blocks[1][1]
    off = blocks[0][4]+25*blocks[1][4]
    lower = arb(CHARACTER)
    if mean < lower:
        return True
    if mean > lower:
        min_u.append(off/(mean-lower))
        min_v.append(-off/(mean-lower))
    u = max(x.lower() for x in min_u)
    v = max(x.lower() for x in min_v)
    if not (u > 0 and v > 0):
        return False
    if u*v > 1:
        return True
    cut(w*(u-1/u)-z+s, arb(1), z-w*u)
    cut(-w*(1/v-v)+z-s, arb(-1), -z-w*v)
    alpha, beta = ((1,25) if chart == 0 else (25,1))
    cut(-alpha*h[1],-beta-alpha*h[2],blocks[0][0]+25*blocks[1][0]+alpha*h[0]-arb(TRACE))
    return not poly or cover_orientations(model,split,poly,blocks,h,chart,rows)


def cover(root):
    model = run_path(str(root / "ch39-global-independent-weak-frame-separation.py"))
    split = run_path(str(root / "ch39-left-right-independent-weak-frame-separation.py"))["split_polygon"]
    model["azimuth"]["compile_cover"](root)
    rows = []
    for j,cap,dist,forms in model["up_rays"](2):
        tightened = [(n,sup*arb(CAPS[1]),a,b,c,n2) for n,sup,a,b,c,n2 in forms]
        rows.append((j,cap*arb(CAPS[0]),dist,tightened))
    old, geom, pencil = (model[k] for k in ("lib","geom","pencil"))
    mix, mass, a, b = (model[k] for k in ("V","mass","A","B"))
    for chart,(base,other) in enumerate(((a,b),(b,a))):
        for idx in range(3):
            stack = [(arb(-1),arb(1),0)]
            count = deepest = leaves = 0
            while stack:
                lo,hi,depth = stack.pop()
                t = lo.union(hi)
                count += 1
                deepest = max(depth,deepest)
                lines = old["lines"](base,other,lo,hi,(idx,),geom)
                raw = pencil["physical_blocks"](lines[idx],chart,t,mix,mass,geom,old,(),complex_cross=True) if idx in lines else None
                if raw:
                    blocks = raw[0]
                    tight = model["centred"](old,geom,mix,mass,a,b,chart,idx,lo,hi) if depth > 3 else None
                    if tight:
                        blocks = [tuple(c[k].intersection(v[k]) if k in (0,1,2,4) else v[k] for k in range(6))
                                  for c,v in zip(blocks,tight,strict=True)]
                    if exclude(model,split,blocks,chart,t,rows):
                        leaves += 1
                        continue
                assert depth < 18, (chart,idx,lo,hi)
                mid = (lo+hi)/2
                stack.extend(((lo,mid,depth+1),(mid,hi,depth+1)))
            print("铅笔图、本征支、铅笔盒、成功叶盒、最大深度",chart,idx,count,leaves,deepest,flush=True)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    certify_observables(root)
    cover(root)
