"""第七至九壳三个弱端点的最少站位支撑。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


SUPPORT = (6,7,8,11,14)
LINES = (9,15,15)
ENERGIES = ((177,180),(894,897),(834,837))


def adj(matrix):
    return matrix.conjugate().transpose()


def take_rows(matrix, indices):
    return acb_mat([[matrix[i,j] for j in range(matrix.ncols())] for i in indices])


def take_columns(matrix, indices):
    return adj(take_rows(adj(matrix),indices))


def leading(matrix, k):
    return acb_mat([[matrix[i,j] for j in range(k)] for i in range(k)]).det().real


def exclude_supports(rows, source, tolerance, supports):
    total = rows*source
    count = [0,0]
    for sites in supports:
        small = take_columns(rows,sites)
        c = take_rows(source,sites)
        g = small*adj(small)
        h = adj(c)*c
        neg_f = small*c-total
        block = acb_mat(5,5)
        for i in range(3):
            for j in range(3):
                block[i,j] = g[i,j]
            for j in range(2):
                block[i,j+3] = neg_f[i,j]
                block[j+3,i] = neg_f[i,j].conjugate()
        for i in range(2):
            for j in range(2):
                block[3+i,3+j] = h[i,j]
        # 若 λmin<-ε，阶末或阶四负主子式通常立即作证。
        for i in range(5):
            block[i,i] += tolerance
        if leading(block,5) < 0 or leading(block,4) < 0:
            count[0] += 1
            continue
        # 若 λmin>ε，全部顺序主子式为正。
        for i in range(5):
            block[i,i] -= 2*tolerance
        if all(leading(block,k)>0 for k in range(1,6)):
            count[1] += 1
            continue
        raise AssertionError(("四站 Gram 无法分离",sites))
    return count


def five_sites(shells, source, levels, native, large):
    rows = acb_mat([[basis[i,j] for j in range(19)] for (basis,_),i in zip(shells,LINES,strict=True)])
    s = take_columns(rows,SUPPORT)
    c = take_rows(source,SUPPORT)
    f = rows*source-s*c
    lift = adj(s)*(s*adj(s)).inv()
    def null(seed):
        return seed-lift*(s*seed)
    def norm(v):
        n = (adj(v)*v)[0,0].real.sqrt()
        assert n>0
        return v/n
    u = norm(null(acb_mat([[acb(j+1,2-j)] for j in range(5)])))
    v = null(acb_mat([[acb((j*j+4)%7,2*j-1)] for j in range(5)]))
    v = norm(v-u*(adj(u)*v))
    null_basis = acb_mat([[u[i,0],v[i,0]] for i in range(5)])
    floor,heavy,top = levels
    mp,mz = arb("1.22089e19"),arb("91.1876")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    finite = native["weak_series"](56)
    target = acb_mat([[((arb(a).union(arb(b))/mp-floor)/(heavy-floor)).sqrt(),0] for a,b in ENERGIES])
    particular = lift*(target-f)
    h = adj(c)*c-adj(particular)*particular
    delta = h.det().real.sqrt()
    assert h[0,0].real>0 and delta>0
    # 二阶正定 Hermitian 矩阵的正平方根。
    eye = acb_mat([[1,0],[0,1]])
    root = (h+eye*delta)/(h[0,0].real+h[1,1].real+2*delta).sqrt()
    x = particular+null_basis*root
    out = acb_mat(source)
    for i,station in enumerate(SUPPORT):
        for j in range(2):
            out[station,j] = x[i,j]
    gram_error = adj(x)*x-adj(c)*c
    assert all(z.contains(0) for z in gram_error.entries())
    assert all(z.contains(0) for z in (rows*out-target).entries())
    for (basis,weights),i,ends in zip(shells,LINES,ENERGIES,strict=True):
        proj = basis*out
        amp = min((floor+(heavy-floor)*abs(proj[j,0])**2+(top-floor)*abs(proj[j,1])**2).lower()
                  for j in range(19) if j!=i)
        assert amp**2>p0
        tail = (large(amp**2,p0,p1)-ref).abs_upper()*sum(abs(w).upper() for j,w in enumerate(weights) if j!=i)
        face = [ref+weights[i]*(finite((arb(q)/mz)**2)-ref)+arb(0,tail) for q in ends]
        print("五站低线、其他振幅下界、面值",i+1,*[w.str(13) for w in (amp,*face)],flush=True)
        assert face[0] > 0 > face[1]


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    old = lib("three-shell-rank-one-obstruction")
    old["experiment"](root,fifth,columns,levels)
    shell = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    shell_data = []
    for n in range(2_306_888_100_092,2_306_888_100_095):
        a,b = lib("adjacent-pvm-transport")["spectrum"](shell,6,n)
        _,nodes,pivots,_ = next(v for v in native["shell_weights"](root,n) if v[0]==6)
        assert all(x.overlaps(y) for x,y in zip(a,nodes,strict=True))
        shell_data.append((adj(b),[w.real for w in pivots[2]]))
    ctx.prec = 192
    source = fifth*columns
    five_sites(shell_data,source,levels,native,lib("full-neutral-free-cutoff")["weak_column"])
    coords = shell["geo"]["sites"]
    distances = [sum((a-b)**2 for a,b in zip(coords[i],coords[j],strict=True))
                 for i,j in ((6,11),(6,14),(7,14),(8,14))]
    assert distances==[12,68,68,68]
    cap = sum((arb(q)/arb("1.22089e19")-levels[0])/(levels[1]-levels[0]) for q in old["THRESHOLDS"]).sqrt()
    tolerance = arb("1e-6")
    assert cap < tolerance
    edges = [(i,j) for i,j,d,_ in shell["geo"]["edges"] if d<68]
    remaining = set(range(19))
    sizes = []
    while remaining:
        component = {min(remaining)}
        for _ in range(18):
            for a,b in edges:
                if a in component or b in component:
                    component.update((a,b))
        sizes.append(len(component))
        remaining -= component
    assert sorted(sizes)==[1,2,2,2,4,4,4]
    total = 0
    for ids in product(*old["CANDIDATES"]):
        rows = acb_mat([[basis[i,j] for j in range(19)] for (basis,_),i in zip(shell_data,ids,strict=True)])
        counts = exclude_supports(rows,source,tolerance,combinations(range(19),4))
        print("三线及四站负、正移位证书数",tuple(i+1 for i in ids),counts,flush=True)
        assert sum(counts)==3876
        total += sum(counts)
    assert total==69768


def verify(root):
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    return experiment(root,*source["verify"](root,return_operator="spectral"))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
