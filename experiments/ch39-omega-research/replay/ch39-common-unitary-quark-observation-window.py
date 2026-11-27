"""十维输入盒的共同酉输运全尺度排除，用参考铅笔本征基中的 Schur 补包围本征线。"""
from pathlib import Path
from runpy import run_path
from sys import argv

from flint import acb, acb_mat, arb


def line(model, chart, lo, hi, branch):
    adj = model["geom"]["adj"]
    a, b = model["A"], model["B"]
    base, other = (a,b) if chart == 0 else (b,a)
    mid = (lo+hi)/2
    center = base+other*mid
    h = acb_mat([[center[i,j].mid() if i <= j else center[j,i].mid().conjugate()
                  for j in range(3)] for i in range(3)])
    for i in range(3):
        h[i,i] = h[i,i].real
    values, vectors = h.eig(right=True, algorithm="rump")
    ids = sorted(range(3), key=lambda j: values[j].real.mid())
    if not all(values[i].real < values[j].real for i,j in zip(ids,ids[1:])):
        return None
    norms = [sum(abs(vectors[k,j])**2 for k in range(3)).sqrt() for j in ids]
    q = acb_mat([[vectors[i,j]/n for j,n in zip(ids,norms,strict=True)] for i in range(3)])
    masses = model["mass"]
    d = [model["geom"]["diag"]([x/row[-1] for x in row[:2]] + [arb(1)]) for row in masses]
    vq = model["V"]*q
    a, b = adj(vq)*d[0]*vq, adj(q)*d[1]*q
    first, second = (a,b) if chart == 0 else (b,a)
    h = first+second*mid+second*arb(0,((hi-lo)/2).upper())
    discs = [h[i,i].real+arb(0,sum(abs(h[i,j]) for j in range(3) if j != i).upper()) for i in range(3)]
    rest = [j for j in range(3) if j != branch]
    if any(not (discs[j]<discs[branch] if j<branch else discs[j]>discs[branch]) for j in rest):
        return None
    interval = discs[branch]
    for _ in range(3):
        minor = acb_mat([[h[i,j]-(interval if i == j else 0) for j in rest] for i in rest])
        try:
            v = minor.solve(acb_mat([[-h[i,branch]] for i in rest]))
        except ZeroDivisionError:
            return None
        x = dict(zip(rest,v.entries()))
        interval = interval.intersection((h[branch,branch]+sum(h[branch,i]*x[i] for i in rest)).real)
    v = acb_mat([[x.get(i,acb(1))] for i in range(3)])
    return q*v*(1/sum(abs(z)**2 for z in v.entries()).sqrt())


def certify_branch(model, rows, split, eliminate, chart, branch):
    geom, lib, pencil = (model[x] for x in ("geom","lib","pencil"))
    a,b,mix,mass = (model[x] for x in ("A","B","V","mass"))
    stack = [(arb(-1),arb(1),0)]
    count = deepest = 0
    while stack:
        lo,hi,depth = stack.pop()
        t = lo.union(hi)
        count += 1
        deepest = max(deepest,depth)
        y = line(model,chart,lo,hi,branch)
        raw = pencil["physical_blocks"](y,chart,t,mix,mass,geom,lib,(),complex_cross=True) if y is not None else None
        if raw:
            blocks = raw[0]
            tight = None
            if depth > 3:
                try:
                    tight = model["centred"](lib,geom,mix,mass,a,b,chart,branch,lo,hi)
                except (ValueError,AssertionError,ZeroDivisionError):
                    tight = None
            if tight:
                profiles = []
                for gauge in (tight,blocks):
                    profiles.append([tuple(c[k].intersection(v[k]) if k in (0,1,2,4) else j[k] for k in range(6))
                                     for c,v,j in zip(blocks,tight,gauge,strict=True)])
            else:
                profiles = [blocks]
            if any(eliminate(model,split,q,chart,t,rows) for q in profiles):
                continue
        if depth == 24:
            raise ArithmeticError(("当前输入盒未获全域覆盖",chart,branch,lo,hi))
        mid = (lo+hi)/2
        stack.extend(((lo,mid,depth+1),(mid,hi,depth+1)))
    print("铅笔图、本征支、区间、成功叶盒、深度",chart,branch,count,(count+1)//2,deepest,flush=True)


def verify(root, fraction):
    base = run_path(str(root / "ch39-all-spaces-common-unitary-quark-right.py"))
    base["certify_observables"](root)
    # 中心回放逐块核验了比原门槛高 0.000319·rank(L)·h_u 的谱界。
    # 这三项只由固定源谱、目的磁拉普拉斯决定，可在整个观测盒中复用。
    stronger_trace = ".012504"
    stronger_character = ".001052"
    reserve = arb(".000319")*arb(base["CAPS"][0])
    assert arb(base["TRACE"])+2*reserve > arb(stronger_trace)
    assert arb(base["CHARACTER"])+reserve > arb(stronger_character)
    eliminate = base["exclude"]
    eliminate.__globals__.update(TRACE=stronger_trace, CHARACTER=stronger_character)
    window = run_path(str(root / "ch39-five-shell-flavor-robustness.py"))
    model = window["input_model"](root,arb(fraction))
    split = run_path(str(root / "ch39-left-right-independent-weak-frame-separation.py"))["split_polygon"]
    model["azimuth"]["compile_cover"](root)
    rows = []
    for j,cap,dist,forms in model["up_rays"](2):
        clipped = [(n,s*arb(base["CAPS"][1]),a,b,c,n2) for n,s,a,b,c,n2 in forms]
        rows.append((j,cap*arb(base["CAPS"][0]),dist,clipped))
    for chart in range(2):
        for branch in range(3):
            certify_branch(model,rows,split,eliminate,chart,branch)
    print("共同酉全空间排除的十维联合半宽倍数",fraction,flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent, argv[1] if len(argv)>1 else ".007")
