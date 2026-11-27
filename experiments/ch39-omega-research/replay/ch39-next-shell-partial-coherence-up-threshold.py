"""部分保留 N+2 上型相干时，N+3 任意输运恢复轻上质量的临界比例。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(m):
    return m.conjugate().transpose()


def source_data(root):
    get = lambda name: run_path(str(root/name))
    rank = get("ch39-rank-two-generator-up-quark.py")
    full = get("ch39-one-generator-gravity-full-ckm.py")
    center,radius = full["verify"](root)
    n = rank["N"]
    _,labels,_ = get("ch39-common-embedding-central-right-ckm.py")["verify_nodes"](
        root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    shell = get("ch39-three-ray-chiral-boundary.py")
    spec = get("ch39-adjacent-pvm-transport.py")["spectrum"]
    source,middle,target = [spec(shell,-4,node)[1] for node in (n+1,n+2,n+3)]
    ctx.prec = 384
    eye = rank["eye"](19)
    x0 = get("ch39-one-generator-gravity-up-strange-bottom.py")["cayley"](
        root,rank["eye"],"ch39-full-ckm-gravity-cayley.tsv")
    hc = acb_mat(19,19)
    directions = full["DIRECTIONS"]
    assert len(set(directions))==len(directions)==len(center)
    for step,(i,j,axis) in zip(center,directions,strict=True):
        z = acb(0,step) if axis=="i" else acb(step)
        hc[i,j] += z
        if i!=j:
            hc[j,i] -= z.conjugate()
    x = (2*(eye-hc).inv()-eye)*x0
    radii = [abs((z+arb(0,radius.upper()))-z).upper() for z in center]
    dx = 2*sum(r*r*(1 if i==j else 2) for r,(i,j,_) in zip(radii,directions,strict=True)).sqrt()
    original = labels[n+1,-4]
    assert all(w>0 for w in original) and sum(original)/19>arb(".06")
    q = adj(middle)*adj(x)*adj(x)*source
    A = acb_mat([[q[i,j]*original[j] for j in range(19)] for i in range(19)])*adj(q)
    diagonal = acb_mat([[A[i,i].real if i==j else 0 for j in range(19)] for i in range(19)])
    amplitude_error = 2*dx*max(w.upper() for w in original).sqrt()
    assert amplitude_error<arb("2.271e-8")
    return A,diagonal,target,amplitude_error


def target_data(target):
    root3 = arb(3).sqrt()
    rays = [(target[14,j],sum(target[i,j] for i in (8,7,6))/root3) for j in range(19)]
    keep = [j for j,(a,b) in enumerate(rays) if abs(a)**2+abs(b)**2>0]
    assert keep==[0,3,6,9,11,15,17]
    assert all(a.contains(0) and b.contains(0) for j,(a,b) in enumerate(rays) if j not in keep)
    a,b = zip(*(rays[j] for j in keep),strict=True)
    entries = ([abs(z)**2 for z in a],[abs(z)**2 for z in b],
               [v*w.conjugate() for v,w in zip(a,b,strict=True)])
    chars = []
    omega = acb(-1,root3)/2
    for o in (omega,omega.conjugate()):
        strength = [abs(target[8,j]+o*target[7,j]+o**2*target[6,j])**2/3 for j in range(19)]
        assert sum(strength).contains(1) and sum(p>0 for p in strength)==6
        chars.append(sorted((max(arb(0),p.lower()) for p in strength),reverse=True))
    return entries,chars


def verify(root):
    A,diagonal,target,error = source_data(root)
    entries,chars = target_data(target)
    small_eigen = run_path(str(root/"ch39-next-shell-arbitrary-up-transport-floor.py"))["small_eigen"]
    def minimum(r_text):
        r = arb(r_text)
        op = (1-r)*diagonal+r*A
        levels = op.eig(algorithm="rump")
        assert all(z.imag.contains(0) and z.real>0 for z in levels)
        roots = sorted((z.real.sqrt() for z in levels),key=lambda z:z.mid())
        assert all(z>w for z,w in zip(roots[1:],roots)) and roots[0]>error
        low = [(z-error).lower()**2 for z in roots]
        high = [(z+error).upper()**2 for z in roots]
        char_floor = min(sum(w*p for w,p in zip(low,c,strict=True)) for c in chars)
        lower = min(small_eigen(perm,entries).lower() for perm in permutations(low[:7]))
        upper = min(small_eigen(perm,entries).upper() for perm in permutations(high[:7]))
        assert upper<char_floor
        return lower.union(upper)
    # 每个临界的两侧均用全部十九个本征值及七条平凡目的线的重排来校验。
    for r1,r2,m in ((".52080",".52090",".00000629"),
                    (".56385",".56395",".00000615"),
                    (".60385",".60395",".00000601")):
        earlier,later = minimum(r1),minimum(r2)
        assert earlier>arb(m)>later,(r1,r2,earlier,later)
        print("保留比例、左侧最优质量下界、右侧最优质量上界、判定值",
              r1,r2,earlier.lower().str(13),later.upper().str(13),m,flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
