"""冻结完整见证的 N+2 输出谱后，求 N+3 所有上型酉输运的四站轻质量最优值。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(m):
    return m.conjugate().transpose()


def small_eigen(weights, entries):
    a,d,z = [sum(w*x for w,x in zip(weights,col,strict=True)) for col in entries]
    return (a+d-((a-d)**2+4*abs(z)**2).sqrt())/2


def verify(root):
    get = lambda name: run_path(str(root/name))
    rank = get("ch39-rank-two-generator-up-quark.py")
    full = get("ch39-one-generator-gravity-full-ckm.py")
    center,radius = full["verify"](root)
    n = rank["N"]
    _, labels, _ = get("ch39-common-embedding-central-right-ckm.py")["verify_nodes"](
        root, ((n+1,".0294"),), ".1", (), ((6,7,10,13),(6,7,10,12)))
    shell = get("ch39-three-ray-chiral-boundary.py")
    spec = get("ch39-adjacent-pvm-transport.py")["spectrum"]
    source, trained, target = [spec(shell,-4,j)[1] for j in (n+1,n+2,n+3)]
    ctx.prec = 384
    eye = rank["eye"](19)
    x0 = get("ch39-one-generator-gravity-up-strange-bottom.py")["cayley"](
        root,rank["eye"],"ch39-full-ckm-gravity-cayley.tsv")
    h = acb_mat(19,19)
    directions = full["DIRECTIONS"]
    for step,(i,j,part) in zip(center,directions,strict=True):
        a = acb(0,step) if part=="i" else acb(step)
        h[i,j] += a
        if i!=j:
            h[j,i] -= a.conjugate()
    assert len(set(directions))==len(directions)==18
    x = (2*(eye-h).inv()-eye)*x0
    radii = [abs((z+arb(0,radius.upper()))-z).upper() for z in center]
    dx = 2*sum(r*r*(1 if i==j else 2) for r,(i,j,_) in zip(radii,directions,strict=True)).sqrt()
    original = labels[n+1,-4]
    weights = rank["read"](source,trained,original,adj(x)*adj(x))
    amplitude_error = 2*dx*max(w.upper() for w in original).sqrt()
    boxes = [(max(arb(0),w.sqrt()-amplitude_error)**2).union((w.sqrt()+amplitude_error)**2)
             for w in weights]
    low = sorted(w.lower() for w in boxes)
    high = sorted(w.upper() for w in boxes)
    assert (shell["soft"]["pole"], *shell["boundary_sites"]()) == (14,8,7,6)
    root3 = arb(3).sqrt()
    r = [[target[14,j], sum(target[i,j] for i in (8,7,6))/root3] for j in range(19)]
    mids = [i for i,(a,b) in enumerate(r) if abs(a)**2+abs(b)**2>0]
    assert len(mids)==7 and all(a.contains(0) and b.contains(0) for i,(a,b) in enumerate(r) if i not in mids)
    omega = acb(-1,root3)/2
    bounds = []
    for o in (omega,omega.conjugate()):
        strength = [abs(target[8,j]+o*target[7,j]+o**2*target[6,j])**2/3 for j in range(19)]
        assert sum(strength).contains(1) and sum(p>0 for p in strength)==6
        mins = sorted((max(arb(0),p.lower()) for p in strength),reverse=True)
        maxs = sorted((p.upper() for p in strength),reverse=True)
        floor = sum(w*v for w,v in zip(low,mins,strict=True))
        ceiling = sum(w*v for w,v in zip(high,maxs,strict=True))
        bounds.append(floor.union(ceiling))
    rays = [r[i] for i in mids]
    entries = ([abs(a)**2 for a,b in rays],[abs(b)**2 for a,b in rays],[a*b.conjugate() for a,b in rays])
    lower = min(small_eigen(perm,entries).lower() for perm in permutations(low[:7]))
    upper = min(small_eigen(perm,entries).upper() for perm in permutations(high[:7]))
    best = lower.union(upper)
    assert min(v.lower() for v in bounds)>arb(".0000260") and best < min(v.lower() for v in bounds)
    assert arb(".000007178") < best < arb(".000007181")
    threshold = arb(".00000629")
    assert best-threshold > arb(".000000888")
    print("两个非平凡字符各自可达的最小值包围",[(v.lower().str(15),v.upper().str(15)) for v in bounds])
    print("平凡块七个目的层（从一开始编号）", [j+1 for j in mids])
    print("第 N+3 壳的全局最小上型质量包围、超过参考窗上端的量",best.lower().str(15),best.upper().str(15),
          (best-threshold).lower().str(15))
    return best,bounds


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
