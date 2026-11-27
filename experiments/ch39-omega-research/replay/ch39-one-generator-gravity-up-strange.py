"""同一十九维发生元实现荷六三读数、上型三质量中心及奇夸克参考窗。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(x):
    return x.conjugate().transpose()


def frob(x):
    return sum(abs(z).upper()**2 for z in x.entries()).sqrt().upper()


def cayley(root, eye):
    data = (root/"ch39-gravity-strange-seven-cayley.tsv").read_text().splitlines()
    h = acb_mat(19,19)
    seen = set()
    for line in data:
        i,j,a,b = map(int,line.split())
        assert 0<=i<=j<19 and (i,j) not in seen and (i!=j or b==0)
        seen.add((i,j))
        z = acb(arb(a)/10**13,arb(b)/10**13)
        h[i,j],h[j,i] = z,z.conjugate()
    assert len(seen)==190
    unit = eye(19)
    return 2*(unit-h*acb(0,1)).inv()-unit


def bracket(value, norm):
    root = value.sqrt()
    assert root > norm
    return ((root-norm)**2).union((root+norm)**2)


def four_block(root, spectrum, shells, n, charge, labels, x):
    rank = run_path(str(root/"ch39-rank-two-generator-up-quark.py"))
    old,target = [spectrum(shells,charge,k)[1] for k in (n,n+1)]
    weights = rank["read"](old,target,labels[n,charge],x)
    four = acb_mat([[target[k,j] for j in range(19)] for k in (14,8,7,6)])
    return acb_mat([[four[i,j]*weights[j] for j in range(19)] for i in range(4)])*adj(four)


def verify(root):
    rank = run_path(str(root/"ch39-rank-two-generator-up-quark.py"))
    spec = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    shells = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    joint = run_path(str(root/"ch39-one-generator-up-and-gravity.py"))
    high = run_path(str(root/"ch39-refine-two-gravity-images.py"))["high_images"]
    s,q,y,z = high(root)
    n = rank["N"]+1
    train = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    _,labels,_ = train(root,((n,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    x = cayley(root,rank["eye"])
    g = x*x*x
    cols = joint["side_by_side"]
    mismatch = frob(g*cols([s,q])-cols([y,z]))
    # 配准两副等距列的极分解给 ||G-G0||_F ≤ 4||[y,z]-G0[s,q]||_F。
    dg = 4*mismatch
    roots = x.eig(algorithm="rump")
    gap = min(abs(u**3-v**3).lower() for j,u in enumerate(roots) for v in roots[j+1:])
    assert gap>2*dg
    coef = min(abs(u**2+u*v+v**2).lower() for u in roots for v in roots)
    # 在 Frobenius 球上对立方映射作压缩，线性逆的范数恰为 1/coef。
    dx = 2*dg/coef
    assert (dg+3*dx**2+dx**3)/coef<dx and (6*dx+3*dx**2)/coef<1
    print("高源像偏差、立方逆系数、发生元误差上界",mismatch.str(10),coef.str(10),dx.str(10),flush=True)

    up = four_block(root,spec,shells,n,-4,labels,adj(x)*adj(x))
    down = four_block(root,spec,shells,n,2,labels,x)
    eps = [max(v.upper() for v in labels[n,c]).sqrt()*err for c,err in ((-4,2*dx),(2,dx))]
    sqrt3 = arb(3).sqrt()
    omega = acb(-arb(1)/2,sqrt3/2)
    chars = [acb_mat([[0],[1],[o],[o**2]])*(1/sqrt3) for o in (omega,omega.conjugate())]
    trivial = acb_mat([[1,0],[0,1/sqrt3],[0,1/sqrt3],[0,1/sqrt3]])
    ub = adj(trivial)*up*trivial
    eig,vec = ub.eig(right=True,algorithm="rump")
    ix = sorted(range(2),key=lambda i:eig[i].real.mid())
    u0,u3 = [eig[i].real for i in ix]
    uc = [(adj(v)*up*v)[0,0].real for v in chars]
    which = sorted(range(2),key=lambda i:uc[i].mid())
    ref = [u0,uc[which[0]],u3,uc[which[1]]]
    eigen = [bracket(t,eps[0]) for t in ref]
    mass = list(map(arb,rank["CENTERS"]))
    assert all(eigen[i]<v<eigen[i+1] for i,v in enumerate(mass))
    probs = [prod(a-v for v in mass)/prod(a-b for j,b in enumerate(eigen) if j!=i) for i,a in enumerate(eigen)]
    assert all(v>0 for v in probs)
    print("上型四站谱的严格包围",[v.str(11) for v in eigen],flush=True)

    # 用实际高平凡本征矢 v 把两字符线精确提升到三平面 n⊥。
    center_vec = acb_mat([[vec[i,ix[1]]] for i in range(2)])
    center_vec *= 1/(adj(center_vec)*center_vec)[0,0].real.sqrt()
    db = adj(trivial)*down*trivial
    mid_value = (adj(center_vec)*db*center_vec)[0,0].real
    op_err = 4*max(v.upper() for v in labels[n,-4])*dx
    sin_bound = op_err/(u3-u0-op_err)
    assert sin_bound>0 and sin_bound<arb(".001")
    max_d = max(v.upper() for v in labels[n,2])
    beta = mid_value+2*max_d*sin_bound+2*max_d*dx
    dc = [bracket((adj(chars[k])*down*chars[k])[0,0].real,eps[1]) for k in which]
    ds = [bracket(v.real,eps[1]) for v in db.eig(algorithm="rump")]
    low,high_edge = arb(".000265"),arb(".000271")
    assert all(v>low for v in [*dc,*ds]) and all(v<high_edge for v in dc)
    # 对任何新三平面和任何整体相位，最轻质量的平方根均服从同一扰动界。
    separation = (low.sqrt()-arb(".0000135").sqrt())/max_d.sqrt()
    assert separation>arb(".02325")
    print("命中下夸克中心所需的发生元射影距离下界",separation.lower().str(12),flush=True)
    loss = (beta-high_edge)*sum(probs[j]/(high_edge-v) for j,v in zip((1,3),dc,strict=True))
    assert loss < probs[2], (loss,probs[2],beta,dc)
    print("下型两字符包围、高平凡上型轴的下型均值上界、二平面负定余量",
          [v.str(11) for v in dc],beta.upper().str(11),(probs[2]-loss).lower().str(11),flush=True)
    # 低平凡上型本征线投回三平面。原射线与归一化投影的正弦距离恰为 sqrt(p_1)。
    low_vec = acb_mat([[vec[i, ix[0]]] for i in range(2)])
    low_vec *= 1/(adj(low_vec)*low_vec)[0,0].real.sqrt()
    low_expectation = (adj(low_vec)*db*low_vec)[0,0].real
    bottom = low_expectation - 2*max_d*(sin_bound+dx) - max_d*probs[0].sqrt()
    assert bottom > arb(".0975")
    print("所有上型中心三平面的底夸克严格下界",bottom.lower().str(12),flush=True)
    print("三个上型中心与奇夸克开区间",*[v.str(8) for v in mass],low.str(8),high_edge.str(8),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
