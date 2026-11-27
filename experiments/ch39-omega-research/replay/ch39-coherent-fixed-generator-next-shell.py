"""右壳见证保留中间壳相干时，固定生成元、固定三平面的下一壳读数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def verify(root):
    def load(name):
        return run_path(str(root/name))
    rank = load("ch39-rank-two-generator-up-quark.py")
    full = load("ch39-one-generator-gravity-full-ckm.py")
    holdout = load("ch39-full-ckm-generator-next-shell.py")
    embeddings = load("ch39-common-embedding-central-right-ckm.py")
    geom = load("ch39-independent-unitary-quark-right-holdout.py")["GEOMETRY"]
    shell = load("ch39-three-ray-chiral-boundary.py")
    spec = load("ch39-adjacent-pvm-transport.py")["spectrum"]
    win = load("ch39-three-shell-outer-holdout.py")["WINDOWS"]
    center,radius = full["verify"](root)
    n = rank["N"]
    _,labels,_ = embeddings["verify_nodes"](
        root, ((n+1,".0294"),), ".1", (), ((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    x0 = load("ch39-one-generator-gravity-up-strange-bottom.py")["cayley"](
        root,rank["eye"],"ch39-full-ckm-gravity-cayley.tsv")
    eye = rank["eye"](19)
    h = acb_mat(19,19)
    directions = full["DIRECTIONS"]
    assert len(set(directions))==len(center)==18
    for t,(i,j,part) in zip(center,directions,strict=True):
        z = acb(t) if part=="r" else acb(0,t)
        h[i,j] += z
        if i!=j:
            h[j,i] -= z.conjugate()
    x = (2*(eye-h).inv()-eye)*x0
    radii = [abs((t+arb(0,radius.upper()))-t).upper() for t in center]
    dx = 2*sum(rt**2*(1 if i==j else 2) for rt,(i,j,_) in zip(radii,directions,strict=True)).sqrt()
    assert dx < arb("1.167e-8")
    keep = (shell["soft"]["pole"],*shell["boundary_sites"]())
    assert keep == (14,8,7,6)
    a,b,c = [1/arb(j).sqrt() for j in (3,2,6)]
    basis = acb_mat([[1,0,0,0],[0,a,b,c],[0,a,-b,c],[0,a,0,-2*c]])
    _,_,embedding,rotation = embeddings["common_embedding"](*geom)
    frame = basis*rotation*embedding
    angles, eigenvectors = [], []
    uniform = []
    def edges(values):
        return [(v.lower().str(11),v.upper().str(11)) for v in values]
    for iq,(charge,op,du) in enumerate(((-4,adj(x)*adj(x),2*dx),(2,x,dx))):
        source,middle,target = [spec(shell,charge,j)[1] for j in (n+1,n+2,n+3)]
        w = labels[n+1,charge]
        spread = max(v.upper() for v in w)-min(v.lower() for v in w)
        err = 2*spread*du
        partial = rank["read"](source,middle,w,op)
        at_zero = rank["read"](middle,target,partial,op)
        at_one = rank["read"](source,target,w,op*op)
        sites = acb_mat([[target[i,j] for j in range(19)] for i in keep])
        cf = adj(frame)*sites
        def projected(weights):
            return acb_mat([[cf[i,j]*weights[j] for j in range(19)] for i in range(3)])*adj(cf)
        a0,a1 = projected(at_zero),projected(at_one)
        values,axes = holdout["sorted_eigen"](a1)
        actual = [v+arb(0,err.upper()) for v in values]
        assert all(v<arb(lo) or v>arb(hi) for v,(lo,hi) in zip(actual,win[iq],strict=True))
        # 每个目的层的精确振幅沿 r 仿射，故全十九层的统一端点界也覆盖任意三平面。
        floor = min(z.lower() for cut in (at_zero,at_one) for z in cut)-err
        cap = max(z.upper() for cut in (at_zero,at_one) for z in cut)+err
        assert floor > arb((".00117514",".00095667")[iq])
        if iq==0:
            assert cap < arb(".644666")
            source_cap = arb(".947555")
            assert max(w.upper() for w in labels[n+1,charge]) < source_cap
            target_edge = arb(".8573")
            angle = (target_edge/source_cap).sqrt().asin()-(arb(".644666")/source_cap).sqrt().asin()
            second_only = 2*(angle/4).sin()
            shared = 2*(angle/8).sin()
            assert second_only > arb(".1434485") and shared > arb(".0717704")
            print("为使顶质量入窗，改第二步或两步共改时生成元至少相距",
                  second_only.lower().str(12),shared.lower().str(12),flush=True)
        # 使用完全相干时的本征标架作为固定试探子空间，端点的下界沿仿射段保持。
        top_two = acb_mat([[axes[1][i,0],axes[2][i,0]] for i in range(3)])
        mid_floor = min(min(z.real.lower() for z in (adj(top_two)*mat*top_two).eig(algorithm="rump"))
                        for mat in (a0,a1))-err
        top_floor = min((adj(axes[2])*mat*axes[2])[0,0].real.lower() for mat in (a0,a1))-err
        if iq==0:
            assert mid_floor>arb(".00316")
        else:
            assert top_floor>arb(".01410")
        uniform.append((mid_floor,top_floor))
        gap = [min(abs(v-u).lower() for j,u in enumerate(values) if j!=i) for i,v in enumerate(values)]
        assert min(gap)>2*err
        angles.append([arb(2).sqrt()*err/(g-err) for g in gap])
        eigenvectors.append(axes)
        print("荷、完全相干时三质量区间",charge,edges(actual),flush=True)
        print("荷、全部相干比例下十九层振幅的统一下上界",charge,floor.lower().str(12),cap.upper().str(12),flush=True)
    overlaps = [[(adj(a)*b)[0,0] for b in eigenvectors[1]] for a in eigenvectors[0]]
    errors = [[(a+b).upper() for b in angles[1]] for a in angles[0]]
    absolute = [[abs(z)+arb(0,e) for z,e in zip(row,er,strict=True)] for row,er in zip(overlaps,errors,strict=True)]
    assert all(z>0 for row in absolute for z in row)
    # 三个标准混合角由三个模唯一恢复。
    s13 = absolute[0][2]
    c13 = (1-s13**2).sqrt()
    sines = [absolute[0][1]/c13,absolute[1][2]/c13,s13]
    references=[arb(m)+arb(0,arb(hr).upper()) for m,hr in
                ((".22704",".00082"),(".04275",".00042"),(".00377",".00008"))]
    assert sines[0]>references[0].sin() and sines[1]<references[1].sin() and sines[2]<references[2].sin()
    print("全比例的固定三平面粲、底质量下界",uniform[0][0].lower().str(12),uniform[1][1].lower().str(12))
    print("完全相干时三个 CKM 角的正弦区间",edges(sines))
    indices=((0,1),(1,2),(0,2),(1,1))
    z = [overlaps[i][j]+acb(arb(0,errors[i][j]),arb(0,errors[i][j])) for i,j in indices]
    signed = (z[0]*z[1]*z[2].conjugate()*z[3].conjugate()).imag
    assert arb("5.09e-8") < signed < arb("5.32e-8")
    print("完全相干时九个 CKM 模平方的区间", [edges([z*z for z in row]) for row in absolute],flush=True)
    print("完全相干时带号 J",edges([signed]),flush=True)
    return absolute,signed


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
