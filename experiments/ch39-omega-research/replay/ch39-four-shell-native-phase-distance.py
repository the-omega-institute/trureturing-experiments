"""连续四壳弱零点到任意第五壳原生谱相位的统一算子距离。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    previous = lib("three-shell-rank-one-obstruction")
    previous["experiment"](root,fifth,columns,levels)
    shell = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    candidates = (*previous["CANDIDATES"],(0,15))
    thresholds = (*previous["THRESHOLDS"],3000)
    native_rows = []
    for n in range(2_306_888_100_092,2_306_888_100_096):
        nodes,basis = lib("adjacent-pvm-transport")["spectrum"](shell,6,n)
        native_rows.append(adj(basis))
    _,ordered,piv,_ = next(v for v in native["shell_weights"](root,n) if v[0]==6)
    assert all(a.overlaps(b) for a,b in zip(nodes,ordered,strict=True))
    weights = [v.real for v in piv[2]]
    assert [i for i,w in enumerate(weights) if w>0]==[0,9,12,15]
    floor,heavy,top = levels
    mp,mz = arb("1.22089e19"),arb("91.1876")
    ctx.prec = 256
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    ref = arb("-.074256577194663365545")
    finite = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    lower,upper = arb("-1.4e-11"),arb(".000016")
    join = arb(10)**10
    assert lower < large(join*p0,p0,p1)-ref < upper
    cells = [((arb(3000)/mz)**2,join)]
    pieces = 0
    while cells:
        a,b = cells.pop()
        value = finite(a.union(b))-ref if b<2*a else None
        if value is not None and lower<value<upper:
            pieces += 1
            continue
        assert b-a > arb("1e-8")
        mid = ((a+b)/2).mid()
        cells.extend(((a,mid),(mid,b)))
    no_low = ref+sum(w*(upper if i in candidates[-1] else (arb(".020483945") if w>0 else lower))
                       for i,w in enumerate(weights))
    print("第十壳两候选均达 3000 GeV 时的弱差上界与分段数",no_low.str(14),pieces,flush=True)
    assert no_low<0
    caps = [((arb(q)/mp-floor)/(heavy-floor)).sqrt().upper() for q in thresholds]
    records = {}
    for count in (3,4):
        epsilon = sum(x**2 for x in caps[:count]).sqrt().upper()
        group = {}
        for picks in product(*candidates[:count]):
            rows = acb_mat([[basis[i,j] for j in range(19)] for basis,i in zip(native_rows[:count],picks,strict=True)])
            gram_inv = (rows*adj(rows)).inv()
            # Q 先由高精度目标谱行构造，随后转入第五壳谱基。
            q = adj(rows)*gram_inv*rows
            projected = adj(fifth)*q*fifth
            # 三行/四行的伪逆算子范数不大于其 Frobenius 范数。
            error = epsilon*sum(gram_inv[i,i].real for i in range(count)).sqrt().upper()
            bound = arb(0)
            for j in range(19):
                diagonal = projected[j,j].real
                if diagonal.lower()<=arb("1e-14"):
                    continue
                norm = diagonal.sqrt()
                for h in range(2):
                    main = norm*abs(columns[j,h])
                    other = sum(abs(projected[j,k]*columns[k,h]) for k in range(19) if j!=k)/norm
                    bound = max(bound,(main-other-error).lower())
            group[picks] = bound
        print("候选壳数、组合数、全部下界、指定前三线下界",count,len(group),min(group.values()).str(15),
              min(v for k,v in group.items() if k[:3]==(0,0,15)).str(15),flush=True)
        records[count] = group
    assert min(records[3].values())>arb(".1838")
    assert min(records[4].values())>arb(".4453")
    assert records[3][0,0,15]>arb(".7671")
    assert min(v for k,v in records[4].items() if k[:3]==(0,0,15))>arb(".8581")
    return records


if __name__=="__main__":
    root = Path(__file__).resolve().parent
    original = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*original["verify"](root,return_operator="spectral"))
