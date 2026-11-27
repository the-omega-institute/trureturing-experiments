"""冻结完整 CKM 见证的整个 Cayley 盒，在下一实际零点复算夸克留出读数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def sorted_eigen(matrix):
    values, vecs = matrix.eig(right=True,algorithm="rump")
    order = sorted(range(len(values)),key=lambda i:values[i].real.mid())
    assert all(z.imag.contains(0) for z in values)
    levels = [values[j].real for j in order]
    assert all(l<r for l,r in zip(levels,levels[1:]))
    cols = [acb_mat([[vecs[i,j]] for i in range(matrix.nrows())]) for j in order]
    return levels,[v*(1/(adj(v)*v)[0,0].real.sqrt()) for v in cols]


def verify(root):
    get = lambda name: run_path(str(root/name))
    full = get("ch39-one-generator-gravity-full-ckm.py")
    center,radius = full["verify"](root)
    rank = get("ch39-rank-two-generator-up-quark.py")
    emb = get("ch39-common-embedding-central-right-ckm.py")
    indep = get("ch39-independent-unitary-quark-right-holdout.py")
    shell = get("ch39-three-ray-chiral-boundary.py")
    spectrum = get("ch39-adjacent-pvm-transport.py")["spectrum"]
    windows = get("ch39-three-shell-outer-holdout.py")["WINDOWS"]
    n = rank["N"]
    _,labels,_ = emb["verify_nodes"](root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    eye = rank["eye"](19)
    x0 = get("ch39-one-generator-gravity-up-strange-bottom.py")["cayley"](root,rank["eye"],"ch39-full-ckm-gravity-cayley.tsv")
    hc = acb_mat(19,19)
    dirs = full["DIRECTIONS"]
    assert len(set(dirs))==len(dirs)==len(center)
    for a,(i,j,axis) in zip(center,dirs,strict=True):
        direction = acb(a) if axis=="r" else acb(0,a)
        hc[i,j] += direction
        if i!=j:
            hc[j,i] -= direction.conjugate()
    x = (2*(eye-hc).inv()-eye)*x0
    # 两个 Cayley 逆均为压缩。所列反厄米方向关于实 Frobenius 内积正交。
    radii = [abs((z+arb(0,radius.upper()))-z).upper() for z in center]
    dx = 2*sum(r*r*(1 if i==j else 2) for r,(i,j,_) in zip(radii,dirs,strict=True)).sqrt()
    assert dx<arb("1.167e-8")
    keep = (shell["soft"]["pole"],*shell["boundary_sites"]())
    assert keep==(14,8,7,6)
    s,t,u = [1/arb(j).sqrt() for j in (3,2,6)]
    four_basis = acb_mat([[1,0,0,0],[0,s,t,u],[0,s,-t,u],[0,s,0,-2*u]])
    _,_,intrinsic,rotation = emb["common_embedding"](*indep["GEOMETRY"])
    frame = four_basis*rotation*intrinsic
    print("覆盖十八参数盒的发生元算子距离上界",dx.upper().str(11),flush=True)
    min_angle, vectors, actual = [],[],[]
    for iq,(charge,op,eps) in enumerate(((-4,adj(x)*adj(x),2*dx),(2,x,dx))):
        source,middle,dest = [spectrum(shell,charge,node)[1] for node in (n+1,n+2,n+3)]
        weights = labels[n+1,charge]
        spread = max(y.upper() for y in weights)-min(y.lower() for y in weights)
        first = rank["read"](source,middle,weights,op)
        second = rank["read"](middle,dest,first,op)
        spread_mid = max(y.upper() for y in first)-min(y.lower() for y in first)
        # 每一次去相干和固定等距压缩都不增大自伴算子误差。
        err = (spread+spread_mid)*eps
        wiggle = arb(0,err.upper())
        sites = acb_mat([[dest[i,j] for j in range(19)] for i in keep])
        raw = acb_mat([[sites[i,j]*second[j] for j in range(19)] for i in range(4)])*adj(sites)
        four,_ = sorted_eigen(raw)
        interlace = [v+wiggle for v in four]
        assert all(interlace[j]>arb(hi) or interlace[j+1]<arb(lo)
                   for j,(lo,hi) in enumerate(windows[iq]))
        levels,vec = sorted_eigen(adj(frame)*raw*frame)
        mass = [v+wiggle for v in levels]
        assert all(val<arb(lo) or val>arb(hi) for val,(lo,hi) in zip(mass,windows[iq],strict=True)),mass
        gap = [min(abs(y-z).lower() for j,z in enumerate(levels) if j!=i) for i,y in enumerate(levels)]
        assert all(g>2*err for g in gap)
        # 投影到中心本征线的正交补，用谱隙去除其本征向量方程的残差。
        min_angle.append([arb(2).sqrt()*err/(g-err) for g in gap])
        vectors.append(vec)
        actual.append(mass)
        print("荷、谱误差、四站谱、三质量",charge,err.upper().str(11),
              [(v+wiggle).str(13) for v in four],[v.str(13) for v in mass],flush=True)
        lower = min(y.lower() for y in second)-err
        upper = max(y.upper() for y in second)+err
        print("荷、第二次输出的全十九层振幅统一下界和上界",charge,lower.str(12),upper.str(12),flush=True)
        # 全十九维振幅去相干后夹在其最小和最大本征值之间。
        if charge==-4:
            assert lower>arb(".00117514") and upper<arb(".644666")<arb(windows[iq][2][0])
        else:
            assert lower>arb(".00106117")>arb(windows[iq][1][1]) and upper<arb(".117555")
        assert lower>arb(windows[iq][0][1])
    overlaps = [[(adj(a)*b)[0,0] for b in vectors[1]] for a in vectors[0]]
    errors = [[(a+b).upper() for b in min_angle[1]] for a in min_angle[0]]
    absolute = [[abs(v)+arb(0,e) for v,e in zip(row,ee,strict=True)] for row,ee in zip(overlaps,errors,strict=True)]
    # 可分别调节两副本征矢的相位进行误差比较；带号 J 在这些相位下不变。
    indices = ((0,1),(1,2),(0,2),(1,1))
    values = []
    for i,j in indices:
        band = arb(0,errors[i][j])
        values.append(overlaps[i][j]+acb(band,band))
    signed = (values[0]*values[1]*values[2].conjugate()*values[3].conjugate()).imag
    angles = [arb(m)+arb(0,arb(e).upper()) for m,e in
              ((".22704",".00082"),(".04275",".00042"),(".00377",".00008"))]
    assert absolute[0][1]>angles[0].sin()*angles[2].cos()
    assert absolute[1][2]<angles[1].sin()*angles[2].cos()
    assert absolute[0][2]<angles[2].sin()
    assert abs(signed)<arb("3e-9")
    print("三个 CKM 模和带号 J",*[absolute[i][j].str(12) for i,j in ((0,1),(1,2),(0,2))],signed.str(12),flush=True)
    return actual,absolute,signed


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
