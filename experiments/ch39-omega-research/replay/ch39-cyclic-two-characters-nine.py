"""三循环等变的同谱算子在两个字符块内共同实现连续三壳九式。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat


PARAMETERS = (
    "-0.625", "0.125", "-3", "-1.625",
    "1", "-0.5", "1.25", "0.375",
    "-1.25", "-0.375", "-0.63862992994126", "-0.375",
    "-0.692973662490878", "0.75", "-1.625", "-1.25",
    "-0.5", "-1.25", "-1.5", "0.375",
    "-1", "-1.375", "2", "-1.875",
    "-1.67703043793728", "-0", "0.125", "0.25",
    "-0.625", "0.875", "-1.75", "1.125",
    "-0.126069857915443", "-0.25", "-0.5", "0",
    "0.375", "-0", "0.875", "0.5",
    "0.073978953524057", "-0.125", "-1.375", "-1",
    "-0.338349403353285", "1.25", "-1.875", "2",
    "0.75", "-0.75", "-1", "-1.75",
)
ACTIVE = (10, 12, 24, 32, 40, 44)
LINES = (0,0,15)
ENERGIES = ("522.398135900","534.208610459","835.615125179")


def adj(a):
    return a.conjugate().transpose()


def eye(n):
    return acb_mat([[int(i==j) for j in range(n)] for i in range(n)])


def characters(sites):
    rotation = [sites.index((b,c,a)) for a,b,c in sites]
    orbits = []
    remaining = set(range(19))
    while remaining:
        start = min(remaining)
        orbit = [start]
        j = rotation[start]
        while j!=start:
            orbit.append(j)
            j = rotation[j]
        remaining.difference_update(orbit)
        orbits.append(orbit)
    omega = acb(-arb(1)/2,arb(3).sqrt()/2)
    frames = []
    for k in range(3):
        selected = [o for o in orbits if k==0 or len(o)==3]
        form = acb_mat(19,len(selected))
        for j,o in enumerate(selected):
            for a,i in enumerate(o):
                form[i,j] = omega**(k*a)/arb(len(o)).sqrt()
        frames.append(form)
    assert [f.ncols() for f in frames]==[7,6,6]
    return frames,rotation


def source_qr(block):
    a,b = [acb_mat([[block[i,j]] for i in range(block.nrows())]) for j in range(2)]
    n = (adj(a)*a)[0,0].real.sqrt()
    assert n>0
    u = a/n
    coeff = (adj(u)*b)[0,0]
    reduced = b-u*coeff
    residual = sum(abs(v)**2 for v in reduced.entries())
    assert residual>0
    return n,coeff,residual


def experiment(root, fifth, columns, levels, *, return_blocks=False):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    general = lib("common-three-shell-nine-endpoints")
    shell = lib("three-ray-chiral-boundary")
    forms,rotation = characters(shell["geo"]["sites"])
    geometry = lib("fifth-seventh-shortest-site-gates")["geometry"]
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    shells = []
    for n,idx in zip(range(2_306_888_100_092,2_306_888_100_095),LINES,strict=True):
        _,_,_,basis,_,weights,heat,field = geometry(root,fifth,columns,levels,number=n)
        # 简单谱的本征行有确定字符。平凡字符分量严格非零，故该行属平凡字符。
        row = acb_mat([[basis[j,idx].conjugate() for j in range(19)]])
        assert any(abs(z)>arb(".01") for z in (row*forms[0]).entries())
        shells.append((adj(basis),weights,heat,field))
    # 荷六连接在原生极点标架下与坐标三循环交换；逐边核对原几何、
    # 对实际两个时钟使用同一正向指标（跨时钟边没有旋转反向）。
    edges = {(a,b):(d,z) for a,b,d,z in shell["geo"]["edges"]}
    assert rotation[shell["soft"]["pole"]]==shell["soft"]["pole"]
    orbit = shell["soft"]["star"]["orbit"]
    assert set(orbit)=={rotation[i] for i in orbit}
    for (a,b),(d,z) in edges.items():
        x,y = rotation[a],rotation[b]
        e,w = edges[tuple(sorted((x,y)))]
        assert d==e and (d not in (724,804) or x<y)
        assert (w if x<y else w.conjugate()).overlaps(z)
    rows = acb_mat([[basis[idx,j] for j in range(19)] for (basis,_,_,_),idx in zip(shells,LINES,strict=True)])
    s = rows*forms[0]
    dual = adj(s)*(s*adj(s)).inv()
    projectors = (eye(7)-dual*s,eye(6))
    original = fifth*columns
    blocks = [adj(f)*original for f in forms]
    factors0 = [source_qr(a) for a in blocks[:2]]
    fixed = forms[2]*blocks[2]
    mp,mz = arb("1.22089e19"),arb("91.1876")
    low,heavy,top = levels
    normalize,inner = general["normalize"],general["inner"]
    def frame(point,energies, return_blocks=False):
        b = acb_mat([[((q/mp-low)/(heavy-low)).sqrt(),0] for q in energies])
        particular = dual*b
        # 以源列的正交残量避免在近共线时直接相减两个 Gram 对角量。
        factors = []
        for k,(n,coeff,residual) in enumerate(factors0):
            loss = sum(abs(particular[i,0])**2 for i in range(7)) if k==0 else arb(0)
            sq = n*n-loss
            remainder = residual-abs(coeff)**2*loss/sq
            assert sq>0 and remainder>0
            diag = sq.sqrt()
            factors.append(acb_mat([[diag,coeff*n/diag],[0,remainder.sqrt()]]))
        ws = list(map(arb,PARAMETERS))
        for j,v in zip(ACTIVE,point,strict=True):
            ws[j] = v
        character_frames = []
        total = [fixed+forms[0]*particular]+[acb_mat(19,2) for _ in point]
        base = 0
        for form,proj,factor in zip(forms[:2],projectors,factors,strict=True):
            n = form.ncols()
            raw = acb_mat([[acb(ws[base+2*i+h],ws[base+2*n+2*i+h]) for h in range(2)] for i in range(n)])
            raw = proj*raw
            ds = [acb_mat(n,2) for _ in point]
            for d,j in zip(ds,ACTIVE,strict=True):
                if base<=j<base+4*n:
                    imag,idx = divmod(j-base,2*n)
                    row,h = divmod(idx,2)
                    for i in range(n):
                        d[i,h] = proj[i,row]*acb(1-imag,imag)
            def col(m,j):
                return acb_mat([[m[i,j]] for i in range(n)])
            a,b = col(raw,0),col(raw,1)
            u,du = normalize(a,[col(d,0) for d in ds])
            c = inner(u,b)
            db = [col(d,1)-z*c-u*(inner(z,b)+inner(u,col(d,1))) for d,z in zip(ds,du,strict=True)]
            v,dv = normalize(b-u*c,db)
            for j,(left,right) in enumerate(zip([u]+du,[v]+dv,strict=True)):
                block = acb_mat([[left[i,0],right[i,0]] for i in range(n)])*factor
                if j==0:
                    character_frames.append(block+particular if base==0 else block)
                total[j] += form*block
            base += 4*n
        return [*character_frames,blocks[2]] if return_blocks else total

    def system(point,qs,values):
        mats = frame(point,qs)
        fs,js = [],[]
        for basis,w,heat,field in shells:
            proj,*ds = [basis*m for m in mats]
            f,j = gravity(proj,ds,levels,w,heat,field,values)
            fs.extend(f.tolist())
            js.extend(j.tolist())
        return arb_mat(fs),arb_mat(js)
    center = [arb(PARAMETERS[i]) for i in ACTIVE]
    q = list(map(arb,ENERGIES))
    eps = arb(".0003")
    qs = [x+arb(0,eps) for x in q]
    for _ in range(3):
        f,j = system(center,q,True)
        step = j.mid().solve(f.mid()).entries()
        center = [(v-d).mid() for v,d in zip(center,step,strict=True)]
    f,j = system(center,qs,True)
    pre = j.mid().inv().mid()
    assert not pre.det().contains(0)
    offsets = [v.abs_upper() for v in (pre*f).entries()]
    widths = [2*d+arb(".000002") for d in offsets]
    box = [c+arb(0,r.upper()) for c,r in zip(center,widths,strict=True)]
    ident = arb_mat([[int(i==j) for j in ACTIVE] for i in ACTIVE])
    norms = [arb(0) for _ in ACTIVE]
    for ix in product(range(4),repeat=6):
        sub = [c+(2*i-3)*r/4+arb(0,(r/4).upper()) for c,r,i in zip(center,widths,ix,strict=True)]
        _,whole = system(sub,qs,False)
        for k,row in enumerate((ident-pre*whole).tolist()):
            limit = sum(z.abs_upper()*r for z,r in zip(row,widths,strict=True))/widths[k]
            norms[k] = max(norms[k],limit)
    ratio = max(o/r+g for o,r,g in zip(offsets,widths,norms,strict=True))
    print("六维中心偏移、最大盒半宽、自映比",*[v.str(13) for v in (max(offsets),max(widths),ratio)],flush=True)
    assert ratio<arb(".904")
    unit = frame(box,qs)[0]
    weak = lib("native-rh-ray-pvm")["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    character_values = frame(box,qs,return_blocks=True)
    for (basis,w,_,_),idx,qq in zip(shells,LINES,q,strict=True):
        overlaps = [basis*e for e in forms]
        projections = [m*a for m,a in zip(overlaps,character_values,strict=True)]
        lower_bounds = []
        for i in range(19):
            if i==idx:
                continue
            choices = [k for k,m in enumerate(overlaps) if any(abs(m[i,j])>arb(".01") for j in range(m.ncols()))]
            assert len(choices)==1
            pr = projections[choices[0]]
            cutoff = low+sum(g*max(arb(0),abs(pr[i,j]).lower())**2 for j,g in enumerate((heavy-low,top-low)))
            lower_bounds.append(cutoff.lower())
        amin = min(lower_bounds)
        assert amin > p0.sqrt()
        tail = (large(amin**2,p0,p1)-ref).abs_upper()*sum(abs(c).upper() for i,c in enumerate(w) if i!=idx)
        faces = [ref+w[idx]*(weak(((qq+sign*eps)/mz)**2)-ref)+arb(0,tail) for sign in (-1,1)]
        print("低线、余振幅下界、弱差面",idx+1,*[v.str(12) for v in (amin,*faces)],flush=True)
        assert faces[0]>0>faces[1]
    if return_blocks:
        return forms,character_values,levels
    return unit,levels


def verify(root):
    original = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    return experiment(root,*original["verify"](root,return_operator="spectral"))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
