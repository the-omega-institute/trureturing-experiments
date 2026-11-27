"""同一发生元实现荷六三读数、三个上型中心及奇夸克、底夸克两个窗口。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def cayley(root, eye, filename):
    data = (root/filename).read_text().splitlines()
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


def setup(root):
    prior = run_path(str(root/"ch39-one-generator-gravity-up-strange.py"))
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
    return prior,rank,spec,shells,joint,n,labels,s,q,y,z


def verify_case(root, data, filename, endpoints, expected_signs):
    prior,rank,spec,shells,joint,n,labels,s,q,y,z = data
    adj,frob,bracket,four_block = [prior[k] for k in ("adj","frob","bracket","four_block")]
    x = cayley(root,rank["eye"],filename)
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


    db = adj(trivial)*down*trivial
    up_move = (2*u3.sqrt()+eps[0])*eps[0]
    sin_bound = up_move/(u3-u0-up_move)
    assert sin_bound>0 and sin_bound<arb(".001")
    mu = [bracket(z.real,eps[1]) for z in db.eig(algorithm="rump")]
    down_move = (2*max(z.upper() for z in mu).sqrt()+eps[1])*eps[1]
    diagonal = []
    for j in ix:
        v = acb_mat([[vec[i,j]] for i in range(2)])
        v *= 1/(adj(v)*v)[0,0].real.sqrt()
        central = (adj(v)*db*v)[0,0].real
        move = down_move+2*max(z.upper() for z in mu)*sin_bound
        diagonal.append(central+arb(0,move.upper()))
    c1,c2 = [bracket((adj(chars[k])*down*chars[k])[0,0].real,eps[1]) for k in which]

    p_low,p_1,p_high,p_2 = probs
    b_low,b_high = diagonal
    # 令平凡平面两个法向坐标间的相位使交叉余子式实部为零。
    def characteristic(t):
        pair = (t-c1)*(t-c2)
        trivial_terms = (p_low*(t-b_high)+p_high*(t-b_low))*pair
        char_terms = (p_1*(t-c2)+p_2*(t-c1))*prod(t-z for z in mu)
        return trivial_terms+char_terms

    # 六个质量中心的相位检验：法向相对相位只改变 (t-c1)(t-c2) 的系数。
    desired = tuple(map(arb, (".0000135", ".000268", ".01401")))
    residuals = [characteristic(c)-prod(c-m for m in desired) for c in (c1,c2)]
    nonzero = [abs(r).lower() for r in residuals]
    assert max(nonzero)>0
    print("在两条下型字符处的质量特征多项式差", [r.str(12) for r in residuals],flush=True)

    points = tuple(map(arb,endpoints))
    signs = tuple(characteristic(t) for t in points)
    assert all(value*sign>0 for value,sign in zip(signs,expected_signs,strict=True)),signs
    print("相位选定后三根的包围端点及特征多项式符号",endpoints,expected_signs,flush=True)
    print("各端点符号检验的最小绝对余量",min(abs(v).lower() for v in signs).str(10),flush=True)
    return sorted([c1,c2,*mu],key=lambda t:t.mid()),max(v.upper() for v in labels[n,2])


def verify(root):
    data = setup(root)
    raw, norm = verify_case(root,data,"ch39-gravity-strange-bottom-cayley.tsv",
                            (".0002669",".0002674",".000271",".01392",".01410"),(-1,1,-1,-1,1))
    assert min(z.lower() for z in raw)>arb(".0002669")
    distance = (arb(".0002669").sqrt()-arb(".0000137").sqrt())/norm.sqrt()
    assert distance>arb(".0233")
    print("第一组见证进入最轻下型窗所需的发生元距离",distance.lower().str(12),flush=True)
    raw, norm = verify_case(root,data,"ch39-gravity-light-down-cayley.tsv",
                            (".0000133",".0000137",".00125",".00130",".094",".096"),(-1,1,1,-1,-1,1))
    assert sum(z>arb(".00126") for z in raw)>=3
    distance = (arb(".00126").sqrt()-arb(".000271").sqrt())/norm.sqrt()
    assert distance>arb(".035")
    print("第二组见证进入奇窗所需的发生元距离",distance.lower().str(12),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
