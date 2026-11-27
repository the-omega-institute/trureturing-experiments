"""按预登记的极点闭路相位，共同复算第十一至十三壳的六个引力读数和固定门对照。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


NINTH = 2_306_888_100_094


def choose_loops(shell):
    pole = shell["soft"]["pole"]
    assert pole == 14
    lengths = {(a,b):d for a,b,d,_ in shell["geo"]["edges"]}
    loops = {}
    for i in range(19):
        if i == pole:
            continue
        for count in (1, 2):
            candidates = []
            for tail in permutations([j for j in range(19) if j not in (i,pole)], count):
                loop = (pole, i, *tail, pole)
                ds = []
                balance = {724:0, 804:0}
                for a,b in zip(loop,loop[1:]):
                    key = tuple(sorted((a,b)))
                    if key not in lengths:
                        break
                    d = lengths[key]
                    ds.append(d)
                    if d in balance:
                        balance[d] += 1 if a < b else -1
                if len(ds) == len(loop)-1 and any(balance.values()):
                    candidates.append((sum(ds),loop,tuple(balance.values())))
            if candidates:
                _, route, exponents = min(candidates)
                loops[i] = route, exponents
                break
        assert i in loops
    return pole, loops


def edge(matrix, a, b):
    z = -matrix[a,b]
    assert abs(z) > 0
    return z/abs(z)


def loop_phase(matrix, route):
    z = acb(1)
    for a,b in zip(route,route[1:]):
        z *= edge(matrix,a,b)
    return z


def check_gauge(matrix, frame, pole, loops, seed):
    units = []
    for i in range(19):
        t = acb(0, arb(i+1)*seed/19)
        units.append((1+t)/(1-t))
    gauge = frame["diagonal"](units)
    changed = gauge * matrix * gauge.conjugate().transpose()
    spokes = frame["spoke_frame"](changed,pole)
    expected = frame["diagonal"]([v*units[pole].conjugate() for v in units])
    assert all(a.overlaps(b) for a,b in zip(spokes.entries(),expected.entries(),strict=True))
    for route,_ in loops.values():
        assert loop_phase(changed,route).overlaps(loop_phase(matrix,route))


def cmp_error(candidate, baseline):
    if candidate < baseline:
        return "闭路"
    if candidate > baseline:
        return "固定"
    return "未分辨"


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    shell = lib("three-ray-chiral-boundary")
    ordered = lib("two-shell-ordered-bond-gates")
    tri = lib("three-shell-gravity-radius56")
    frame = lib("pole-gauge-covariant-transport")
    geometry = lib("fifth-seventh-shortest-site-gates")["geometry"]
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    layers = ordered["LAYERS"] + tri["EXTRA_LAYERS"]
    box = tri["experiment"](root,fifth,columns,levels)
    center = [v.mid() for v in box]
    pole, loops = choose_loops(shell)
    print("极点闭路及 724、804 的有向净计数", loops, flush=True)
    assert sum(len(route) == 4 for route,_ in loops.values()) == 9
    assert sum(len(route) == 5 for route,_ in loops.values()) == 9
    for k,pairs in layers:
        for a,b in pairs:
            da,db = loops[a][1],loops[b][1]
            assert tuple(x-y for x,y in zip(da,db,strict=True)) == ((-2,0) if k in (0,2) else (0,0))
    ctx.prec = 320
    h9 = acb.zeta_zero(NINTH).imag
    m9 = shell["laplacian"](6,h9)
    clocks9 = shell["geo"]["phases"](h9)
    assert set(clocks9) == {724,804}
    basis0 = frame["diagonal"]([1]*19)
    assert frame["spoke_frame"](m9,pole).contains(basis0)
    check_gauge(m9,frame,pole,loops,1)
    initial = [acb(1) if i == pole else loop_phase(m9,loops[i][0]) for i in range(19)]
    targets = arb(1), arb("1e-122")
    errors = []
    for n in (NINTH+2,NINTH+3,NINTH+4):
        h = acb.zeta_zero(n).imag
        m = shell["laplacian"](6,h)
        clocks = shell["geo"]["phases"](h)
        assert set(clocks) == {724,804} and frame["spoke_frame"](m,pole).contains(basis0)
        check_gauge(m,frame,pole,loops,n-NINTH)
        for k,pairs in layers:
            for a,b in pairs:
                assert (edge(m,a,b)*edge(m9,a,b).conjugate()).contains(1)
                if k:
                    assert not edge(m,a,b).imag.contains(0)
        r = [acb(1)]*19
        for i,(route,exponents) in loops.items():
            r[i] = loop_phase(m,route)*initial[i].conjugate()
            simplified = acb(1)
            for d,count in zip((724,804),exponents,strict=True):
                simplified *= (clocks[d]*clocks9[d].conjugate())**(6*count)
            assert abs(r[i]).contains(1) and r[i].overlaps(simplified)
        assert any(not (r[a]*r[b].conjugate()).contains(1) for _,pairs in layers for a,b in pairs)
        _,_,_,basis,_,weights,heat,field = geometry(root,fifth,columns,levels,number=n)
        local_errors = []
        for label,diag in (("固定",basis0),("闭路",frame["diagonal"](r))):
            src = diag.conjugate().transpose()*(fifth*columns)
            left = basis.conjugate().transpose()*diag
            def system(point,values):
                mats = ordered["transport"](src,point,parameters=tri["PARAMETERS"],active=tri["ACTIVE"],layers=layers)
                projection,*directions = [left*x for x in mats]
                return gravity(projection,directions,levels,weights,heat,field,values)
            f,_ = system(center,True)
            _,j = system(box,False)
            residuals = [f[i,0]+sum(j[i,k]*(q-center[k]) for k,q in enumerate(box)) for i in range(2)]
            readings = [v+t for v,t in zip(residuals,targets,strict=True)]
            assert all(not v.contains(0) for v in residuals)
            local_errors.append([abs(v) for v in residuals])
            print(n,label,"牛顿、真空",[(v.lower().str(12),v.upper().str(12)) for v in readings],flush=True)
        comparison = [cmp_error(c,b) for b,c in zip(*local_errors,strict=True)]
        assert comparison == ["固定","固定"]
        errors.append(local_errors)
        print(n,"绝对误差较小者",*comparison,flush=True)
    maxima = [[max(s[m][i].lower() for s in errors).union(max(s[m][i].upper() for s in errors))
               for i in range(2)] for m in range(2)]
    for m,row in enumerate(maxima):
        print(("固定","闭路")[m],"最坏绝对误差",[(v.lower().str(12),v.upper().str(12)) for v in row],flush=True)
    print("最坏绝对误差较小者",*[cmp_error(c,b) for b,c in zip(*maxima,strict=True)],flush=True)
    for row,bounds in zip(maxima, (("1.2664","1.2707",".1753",".1759"),
                                  ("4.5746","4.5831","1.4559","1.4588")), strict=True):
        k,w = row
        lo,hi,wl,wh = map(arb,bounds)
        assert lo < k < hi and wl < w < wh
    return errors,maxima


def verify(root):
    witness = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    return experiment(root,*witness["verify"](root,return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
