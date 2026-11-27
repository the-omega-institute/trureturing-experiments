"""对逐点固定一个六维字符块的三壳九式见证读取第十壳端点。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    forms,blocks,levels,mid,wide,box,center = lib("cyclic-two-characters-nine")["experiment"](
        root,fifth,columns,levels,return_chart=True)
    n = 2_306_888_100_095
    _,_,_,basis,_,weights,heat,field = lib("fifth-seventh-shortest-site-gates")["geometry"](
        root,fifth,columns,levels,number=n)
    overlaps = [basis.conjugate().transpose()*f for f in forms]
    pr = [a*b for a,b in zip(overlaps,blocks,strict=True)]
    projection = acb_mat(19,2)
    for i in range(19):
        choices = [k for k,a in enumerate(overlaps) if any(abs(a[i,j])>arb(".01") for j in range(a.ncols()))]
        assert len(choices)==1
        for j in range(2):
            projection[i,j] = pr[choices[0]][i,j]
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    left = basis.conjugate().transpose()
    f,_ = gravity(left*mid,[],levels,weights,heat,field,True)
    _,derivs = gravity(projection,[left*m for m in wide[1:]],levels,weights,heat,field,False)
    for j in range(2):
        f[j,0] += sum(derivs[j,k]*(q-c) for k,(q,c) in enumerate(zip(box,center,strict=True)))

    p0,p1 = [(arb(v)/arb("1.22089e19"))**2 for v in ("91.1876","1.176397e14")]
    low,heavy,top = levels
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    large = lib("full-neutral-free-cutoff")["weak_column"]
    tail = arb(0)
    assert sum(weights).contains(1)
    for i,w in enumerate(weights):
        lo = low+sum(g*max(arb(0),abs(projection[i,j]).lower())**2 for j,g in enumerate((heavy-low,top-low)))
        assert lo>p0.sqrt()
        tail += abs(w).upper()*(large(lo.lower()**2,p0,p1)-ref).abs_upper()
    values = [f[0,0]+1,f[1,0]+arb("1e-122"),ref+arb(0,tail.upper())]
    print("第十壳牛顿、真空、弱差",[(v.lower().str(15),v.upper().str(15)) for v in values],flush=True)
    for value,(lo,hi) in zip(values,(("1.03247","1.04187"),("-.000695","-.0004715"),
                                    ("-.07425657723","-.07425657716")),strict=True):
        assert arb(lo)<value<arb(hi)
    return values


if __name__=="__main__":
    root = Path(__file__).resolve().parent
    original = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*original["verify"](root,return_operator="spectral"))
