"""右壳同一十九维发生元精确配平荷六三读数与上型三中心，并使三个下型质量同时入窗。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def adj(a):
    return a.conjugate().transpose()


def column(a,j):
    return acb_mat([[a[i,j]] for i in range(a.nrows())])


def frob(a):
    return sum(abs(z).upper()**2 for z in a.entries()).sqrt().upper()


def calibrate(root,x0,rank,spec,shell):
    n = rank["N"]+2
    ordered = run_path(str(root/"ch39-adjacent-ordered-charge-six-cutoffs.py"))["verify"](root)[2][1]
    assert all(z is ordered[0] for z in ordered[:17])
    low,h,k = [ordered[j].sqrt() for j in (0,17,18)]
    gravity = run_path(str(root/"ch39-three-index-free-atoms.py"))
    moments,heat,fields = gravity["broken_fields"](root,[n])
    field = fields[0]
    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    _,_,raw,_ = next(row for row in native["shell_weights"](root,n) if row[0]==6)
    c = [z.real for z in raw[0]]
    assert sum(c).contains(1)
    low_weak = native["weak_series"]()
    high_weak = run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["weak_column"]
    ctx.prec = 384
    mz,mp = arb("91.1876"),arb("1.22089e19")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    low_t,high_t = arb("1.30"),arb("1.50")
    ref = arb("-.074256577194663365545")
    def rough(t):
        return c[15]*low_weak(t)+(1-c[15])*ref
    assert rough(low_t)>0>rough(high_t)
    for _ in range(19):
        mid=(low_t+high_t)/2
        if rough(mid)>0:
            low_t=mid
        else:
            assert rough(mid)<0
            high_t=mid
    full_t = low_t.union(high_t)
    cutoff = (p0*full_t).sqrt()
    eta = ((cutoff-low)/(h-low)).sqrt()
    source,target = [spec(shell,6,j)[1] for j in (n-1,n)]
    high = acb_mat([[source[i,17],source[i,18]] for i in range(19)])
    eye = rank["eye"](19)
    g0 = x0*x0*x0
    old = adj(g0)*column(target,15)
    old_high = adj(high)*old
    phase = old_high[0,0]/abs(old_high[0,0])
    residual = old-high*old_high
    v = column(high,0)*(eta*phase)+residual*((1-eta**2)/(adj(residual)*residual)[0,0].real).sqrt()
    delta = old-v
    # K=old*v*-v*old* 以两列的微小差计算。右乘该旋转会使 G*D16=v。
    K = delta*adj(v)-v*adj(delta)
    cosine=(adj(old)*v)[0,0].real
    assert cosine>arb(".99")
    correction = K+K*K*(1/(1+cosine))
    overlaps = adj(target)*g0*(high+correction*high)
    pairs = [(0,17,acb(0,1)),(2,18,acb(1))]
    source_weights = (h-low,k-low)
    def system(point):
        f = [low+sum(w*abs(overlaps[j,l])**2 for l,w in enumerate(source_weights)) for j in range(19)]
        f[15] = cutoff
        jac = arb_mat(2,2)
        for axis,((j,k_,phase_),angle) in enumerate(zip(pairs,point,strict=True)):
            sn,cs = angle.sin_cos()
            left = [cs*overlaps[j,l]+sn*phase_*overlaps[k_,l] for l in range(2)]
            right = [cs*overlaps[k_,l]-sn*phase_.conjugate()*overlaps[j,l] for l in range(2)]
            f[j],f[k_] = [low+sum(w*abs(q)**2 for w,q in zip(source_weights,row,strict=True)) for row in (left,right)]
            df = 2*sum(w*(a.conjugate()*phase_*b).real for w,a,b in zip(source_weights,left,right,strict=True))
            for i,sign in ((j,1),(k_,-1)):
                y = f[i]
                decay = [(-z/y**2).exp() for _,_,z in field]
                jac[0,axis] += c[i]*sign*df*y*sum(a*b for (_,a,_),b in zip(field,decay,strict=True))/arb.pi()
                jac[1,axis] -= c[i]*sign*df*y**3*sum(g*b for (g,_,_),b in zip(field,decay,strict=True))/(2*arb.pi())
        # 第16层的热积分用绝对界，避免在极小的截止上消去大量有效位。
        err = [abs(c[15])*cutoff**2*sum(abs(a) for _,a,_ in field)/(2*arb.pi()),
               abs(c[15])*cutoff**4*sum(abs(g) for g,_,_ in field)/(8*arb.pi())]
        val = [arb(0,b.upper())-goal for b,goal in zip(err,(arb(1),arb("1e-122")),strict=True)]
        for j,y in enumerate(f):
            if j!=15:
                assert y>arb(".0001")
                terms = moments(field,heat,y**2)
                for i in range(2):
                    val[i] += c[j]*terms[i]
        return arb_mat([[v] for v in val]),jac,f
    point = [arb(0),arb(0)]
    for _ in range(2):
        val,jac,_ = system(point)
        point = [(x-y).mid() for x,y in zip(point,jac.mid().solve(val.mid()).entries(),strict=True)]
    radius = arb("3e-10")
    box = [a+arb(0,radius.upper()) for a in point]
    val,jac,_ = system(point)
    pre = jac.mid().inv()
    _,wide,f = system(box)
    defect = arb_mat([[1,0],[0,1]])-pre*wide
    ratio = max((abs(z)/radius+sum(abs(a) for a in row)).upper() for z,row in zip((pre*val).entries(),defect.tolist(),strict=True))
    assert not pre.det().contains(0) and ratio<1, ratio
    other = sum(c[j]*high_weak(y**2,p0,p1) for j,y in enumerate(f) if j!=15)
    ends = [c[15]*low_weak(t)+other for t in (low_t,high_t)]
    assert ends[0]>0>ends[1],ends
    # 两个正交二阶旋转各有 Frobenius 偏差 ≤ sqrt(2)|angle|。
    dg = frob(correction)+(2*sum(abs(z).upper()**2 for z in box)).sqrt()
    roots = x0.eig(algorithm="rump")
    gap = min(abs(a**3-b**3).lower() for i,a in enumerate(roots) for b in roots[i+1:])
    coef = min(abs(a*a+a*b+b*b).lower() for a in roots for b in roots)
    assert gap>2*dg
    dx = 2*dg/coef
    assert (dg+3*dx**2+dx**3)/coef<dx and (6*dx+3*dx**2)/coef<1
    print("两个旋转中心及共同半宽",*[z.str(9) for z in point],radius.str(3),flush=True)
    print("立方谱隙及其导数最小模下界",gap.str(10),coef.str(10),flush=True)
    print("热方程收缩与自映射界、弱差端点",
          ratio.str(8),ends[0].lower().str(8),ends[1].upper().str(8),flush=True)
    print("第16层 GeV、荷六矩阵改变量、发生元改变量",(mz*full_t.sqrt()).str(12),dg.upper().str(9),dx.upper().str(9),flush=True)
    return dx


def verify(root):
    prior = run_path(str(root/"ch39-one-generator-gravity-up-strange.py"))
    eight = run_path(str(root/"ch39-one-generator-gravity-up-strange-bottom.py"))
    rank = run_path(str(root/"ch39-rank-two-generator-up-quark.py"))
    spec = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    shells = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 384
    x = eight["cayley"](root,rank["eye"],"ch39-gravity-six-quark-cayley.tsv")
    dx = calibrate(root,x,rank,spec,shells)
    n = rank["N"]+1
    train = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    _,labels,_ = train(root,((n,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    bracket,four_block = prior["bracket"],prior["four_block"]
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


    endpoints=(".0000133",".0000137",".000265",".000271",".01392",".01410")
    expected_signs=(-1,1,1,-1,-1,1)
    points = tuple(map(arb,endpoints))
    signs = tuple(characteristic(t) for t in points)
    assert all(value*sign>0 for value,sign in zip(signs,expected_signs,strict=True)),signs
    print("相位选定后三根的包围端点及特征多项式符号",endpoints,expected_signs,flush=True)
    print("各端点符号检验的最小绝对余量",min(abs(v).lower() for v in signs).str(10),flush=True)
    # 压缩后的粲方向正比于 (A-m_c I)^(-1)n。所选法向相位也消去
    # 该方向上下型二次型中的交叉项。以第一字符为基准，避免依赖性损失。
    charm = mass[1]
    relative = [probs[i]/probs[1]*((eigen[1]-charm)/(eigen[i]-charm))**2
                if i != 1 else arb(1) for i in range(4)]
    charm_mean = c1 + sum(w*(value-c1) for w,value in
                          zip(relative,(b_low,c1,b_high,c2),strict=True))/sum(relative)
    # E_c <= s_max + (b_max-s_max) |V_cb|^2。
    vcb_squared = (charm_mean-points[3])/(points[5]-points[3])
    assert vcb_squared > arb(".9967")**2
    print("本构造的粲方向下型期望与 |V_cb| 严格下界",
          charm_mean.str(12),vcb_squared.sqrt().lower().str(10),flush=True)
    # 对任意共同三平面，利用谱投影在该字符线上的权重。
    # 将留数公式中的可消因子先约去，因此端点相等时也可连续使用。
    windows = run_path(str(root/"ch39-common-quadratic-flavor-exclusion.py"))["MASS"]
    up_windows = [arb(lo).union(arb(hi)) for lo,hi in windows[0]]
    m = up_windows[1]
    lam = eigen[1]
    shift = m-lam
    weight = (prod(1+shift/(lam-z) for j,z in enumerate(eigen) if j != 1)
              * prod(1-shift/(m-z) for j,z in enumerate(up_windows) if j != 1))
    # 下型在这条字符线上取 c1，整个正交补也是正半定。
    all_planes = (c1*weight-points[3])/(points[5]-points[3])
    assert all_planes > arb(".95")**2
    print("任意共同三平面、六质量全窗的字符权重与 |V_cb| 下界",
          weight.lower().str(10),all_planes.sqrt().lower().str(10),flush=True)
    # 两端分别向内推进，只保留有严格符号的端点。
    for a,b,sgn in zip(points[::2],points[1::2],expected_signs[::2],strict=True):
        l,r = a,b
        for _ in range(20):
            mid = (l+r)/2
            if sgn*characteristic(mid)>0:
                l = mid
            else:
                r = mid
        lower = l
        l,r = lower,b
        for _ in range(20):
            mid = (l+r)/2
            if sgn*characteristic(mid)<0:
                r = mid
            else:
                l = mid
        assert sgn*characteristic(lower)>0>sgn*characteristic(r)
        print("下型单根开区间",lower.lower().str(14),r.upper().str(14),flush=True)

    source,target = [spec(shells,2,node)[1] for node in (n,n+1)]
    low = [j for j,weight in enumerate(labels[n,2]) if weight < arb(".001")]
    assert len(low)==10 and all(weight>arb(".06") for j,weight in enumerate(labels[n,2]) if j not in low)
    overlap = adj(target)*x*source
    smooth = [sum(abs(overlap[j,i])**2 for i in low).sqrt() for j in range(15,19)]
    assert min(smooth)>dx
    floors = [(a-dx)**2 for a in smooth]
    assert floors[2]>arb(".99999") and floors[0]+floors[2]>arb("1.9998")
    print("第16至19层到最低十个源标签的投影概率下界",[z.lower().str(12) for z in floors],flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
