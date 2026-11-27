"""十八条边同长度共门；为三壳弱差剩余的两组三线模式分别构造精确根。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def experiment(root,fifth,columns,levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    ordered = lib("two-shell-ordered-bond-gates")
    layers = ordered["LAYERS"]+lib("three-shell-gravity-radius56")["EXTRA_LAYERS"]
    graph = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    weak = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    shells = []
    for n in range(2_306_888_100_092,2_306_888_100_095):
        nodes,base = lib("adjacent-pvm-transport")["spectrum"](graph,6,n)
        _,other,pivots,_ = next(r for r in native["shell_weights"](root,n) if r[0]==6)
        assert all(a.overlaps(b) for a,b in zip(nodes,other,strict=True))
        shells.append((adj(base),[z.real for z in pivots[2]]))
    ctx.prec = 256
    eye = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    raw = (root/"ch39-shared-four-gates-weak-boxes.tsv").read_text().splitlines()
    assert len(raw) == 10
    ref = arb("-.074256577194663365545")
    p0,p1 = [(arb(q)/arb("1.22089e19"))**2 for q in ("91.1876","1.176397e14")]
    zero_one = arb(0).union(arb(1))
    for off in (0,5):
        indices = tuple(map(int,raw[off].split()))
        assert indices == ((9,15,15) if off==0 else (15,9,15))
        parameters = raw[off+1].split()
        ds = [[arb(int(x))/10**16 for x in row.split()] for row in raw[off+2:off+5]]
        assert len(parameters)==16 and all(len(d)==16 for d in ds)
        matrices = ordered["transport"](eye,list(map(arb,parameters)),parameters=parameters,active=tuple(range(16)),layers=layers)
        base = matrices[0]
        derivatives = [sum((v*m for v,m in zip(row,matrices[1:],strict=True)),acb_mat(19,19)) for row in ds]
        # 每轮由相同的二阶门组成块对角矩阵。一、二阶坐标偏导的算子范数都至多一。
        # 沿任何直线从原点到立方体中一点，六轮乘积的二阶导数范数至多为 speed²。
        speed = sum(sum(abs(row[j]) for row in ds for j in range(4*k,4*k+4)) for k,_ in layers)
        remainder = arb(0,(speed**2/2).upper())
        noise = acb(remainder,remainder)
        for h,((left,weights),k) in enumerate(zip(shells,indices,strict=True)):
            centered = ((left*base)*fifth)*columns
            dg = [((left*g)*fifth)*columns for g in derivatives]
            def amplitude(i,coords):
                values = [centered[i,c]+sum(t*g[i,c] for t,g in zip(coords,dg,strict=True))+noise for c in range(2)]
                lo = levels[0].lower()+sum((levels[c+1]-levels[0]).lower()*v.abs_lower()**2 for c,v in enumerate(values))
                hi = levels[0].upper()+sum((levels[c+1]-levels[0]).upper()*v.abs_upper()**2 for c,v in enumerate(values))
                return lo.lower(),hi.upper()
            rest = [amplitude(i,[zero_one]*3)[0] for i in range(19) if i!=k]
            minimum = min(rest)
            assert minimum > 0 and minimum**2 > p0 and sum(weights).contains(1)
            tail = sum(w*(large(a**2,p0,p1)-ref) for w,a in
                       zip([w for i,w in enumerate(weights) if i!=k],rest,strict=True))
            faces = []
            for side in (0,1):
                coords = [zero_one]*3
                coords[h] = arb(side)
                a,b = amplitude(k,coords)
                part = weak((a**2/p0).lower().union((b**2/p0).upper()))
                faces.append(ref+weights[k]*(part-ref)+tail)
            print("低线组合、壳、最小余线振幅、弱差两面",indices,h+7,minimum.str(10),faces[0].str(12),faces[1].str(12),flush=True)
            assert faces[0] > 0 > faces[1]


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    src = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*src["verify"](root,return_operator="spectral"))
