"""在同一四站三平面和同一引力达标发生元上扫过底夸克的整个边际窗。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


# 在原生十九站坐标中给一个整数候选。投影进自由空间以后只用它的列空间。
SEED = (
    ((4, -21), (-7, -9), (25, -3), (-19, 8), (-14, -6), (-12, 1)),
    ((2, 7), (13, 9), (10, -38), (-24, 5), (-38, 13), (8, -36)),
    ((-25, -1), (-23, -16), (-12, 3), (12, 3), (4, 19), (25, -11)),
    ((1, 2), (-23, 17), (-17, 11), (6, 4), (11, -22), (-25, 33)),
    ((52, -12), (6, -11), (-6, 5), (14, -8), (-8, 16), (25, 22)),
    ((-28, 5), (31, 10), (-12, -8), (16, -14), (0, -19), (-3, -10)),
    ((8, -2), (5, -10), (-10, 1), (-7, 23), (-6, -29), (5, -6)),
    ((-2, 4), (5, -4), (28, 15), (20, 2), (-21, 6), (-13, -8)),
    ((12, 8), (5, -7), (-3, -36), (6, -1), (15, 2), (10, 5)),
    ((9, 5), (23, 30), (-37, 20), (-20, 3), (-16, 21), (18, 11)),
    ((16, -23), (15, 5), (-8, -29), (20, 23), (3, -33), (-6, -23)),
    ((25, -15), (19, 5), (-11, 2), (9, -27), (16, 7), (-1, -24)),
    ((-21, 10), (23, -11), (2, -10), (20, 12), (-22, 9), (6, 0)),
    ((4, -6), (-40, 2), (19, 9), (23, -16), (-23, 8), (-6, -10)),
    ((12, -3), (-4, -9), (-4, -16), (10, -2), (-17, -7), (-3, 0)),
    ((3, -35), (-8, -14), (2, -16), (-10, 15), (-6, 11), (7, 38)),
    ((6, 5), (-3, 7), (1, 18), (-11, 11), (-17, 15), (-4, -21)),
    ((-19, 22), (18, -33), (13, -11), (4, -12), (21, -6), (16, 21)),
    ((-6, 17), (-19, -5), (-3, 9), (-27, 44), (2, 10), (-1, 1)),
)


def adj(a):
    return a.conjugate().transpose()


def col(a,j):
    return acb_mat([[a[i,j]] for i in range(a.nrows())])


def span(matrix):
    return matrix*(adj(matrix)*matrix).inv()*adj(matrix)


def pair(a,b):
    return acb_mat([[a[i,0],b[i,0]] for i in range(a.nrows())])


def positive(a):
    for k in range(1,a.nrows()+1):
        assert acb_mat([[a[i,j] for j in range(k)] for i in range(k)]).det().real > 0


def verify(root):
    joint = run_path(str(root/'ch39-one-generator-up-and-gravity.py'))
    ranktwo = run_path(str(root/'ch39-rank-two-generator-up-quark.py'))
    orbits = run_path(str(root/'ch39-up-quark-arbitrary-site-phase-floor.py'))['ORBITS']
    train = run_path(str(root/'ch39-common-embedding-central-right-ckm.py'))['verify_nodes']
    lib = run_path(str(root/'ch39-adjacent-pvm-transport.py'))
    sh = run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    _,(a,b),gev,(m,h,_)=joint['verify'](root,return_parameters=True)
    assert arb('103.04')<m*arb('1.22089e19')<arb('103.05')
    assert arb('.4295')<h<arb('.4297')
    n = ranktwo['N']+1
    _,labels,_ = train(root,((n,'.0294'),),'.1',(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    P = ranktwo['projector'](orbits)
    I = ranktwo['eye'](19)
    S = lib['spectrum'](sh,6,n)[1]
    T=lib['spectrum'](sh,6,n+1)[1]
    E=acb_mat([[1/arb(len(g)).sqrt() if i in g else 0 for g in orbits] for i in range(19)])
    s,t = (I-P)*col(S,17),col(S,18)
    _,u0,u1,u2,d,rho=joint['mode_parameters'](P,E,S,T)
    # 第16层回到源第18层的实振幅为 e。由公开源区间和 q16 给其外包。
    e=((gev/arb('1.22089e19')-m)/(h-m)).sqrt()
    def operators(a,b):
        def wave(ee,aa):
            v=u0+ee*u1+aa*u2
            return v+d*(rho-(adj(v)*v)[0,0].real).sqrt()
        # 算朱可比的解析界，再以微小参数半宽乘一次。
        vl=u0+e*u1+a*u2
        rt=(rho-(adj(vl)*vl)[0,0].real).sqrt()
        w=wave(e.mid(),a.mid())
        for direction,var in ((u1,e),(u2,a)):
            w+=(direction-d*((adj(vl)*direction)[0,0].real/rt))*(var-var.mid())
        t3,t19=[joint['align'](col(T,j),t) for j in (2,18)]
        z=t3*(1-b*b).sqrt()+t19*b
        def outer(x):
            return x*adj(x)*(1/(adj(x)*x)[0,0].real)
        def perp(x,v):
            return v-x*((adj(x)*v)[0,0]/(adj(x)*x)[0,0].real)
        ws,zs=perp(s,w),perp(t,z)
        core=P+outer(s)+outer(ws)+outer(t)+outer(zs)
        return I-core,I+(acb(0,1)-1)*P-2*(outer(s-w)+outer(t-z)),w,z
    q=-4
    before = lib['spectrum'](sh,q,n)[1]
    after = lib['spectrum'](sh,q,n+1)[1]
    sites=[sh['soft']['pole'],*sh['boundary_sites']()]
    def block(weights,basis):
        c=acb_mat([[basis[k,j] for j in range(19)] for k in sites])
        return acb_mat([[c[k,j]*weights[j] for j in range(19)] for k in range(4)])*adj(c)
    up = ranktwo['read'](before,after,labels[n,q],I-2*P)
    raw_up=block(up,after)
    # 三个上型中心固定任一候选三平面的法向在两字符线上的概率。
    mu = tuple(map(arb,ranktwo['CENTERS']))
    eigen=sorted([v.real for v in raw_up.eig(algorithm='rump')],key=lambda v:v.mid())
    probs=[prod(l-m for m in mu)/prod(l-r for j,r in enumerate(eigen) if j!=i)
           for i,l in enumerate(eigen)]
    q=2
    before = lib['spectrum'](sh,q,n)[1]
    after = lib['spectrum'](sh,q,n+1)[1]
    seed = acb_mat([[acb(x,y) for x,y in row] for row in SEED])
    weights=labels[n,q]
    def norm(v):
        return sum(abs(x)**2 for x in v.entries()).sqrt()
    free0,old0,w0,z0=operators(a.mid(),b.mid())
    free,old,w,z=operators(a,b)
    def norm_upper(v):
        return sum(abs(x).upper()**2 for x in v.entries()).sqrt().upper()
    dw,dz=norm_upper(w-w0),norm_upper(z-z0)
    def orthogonal(x,v):
        return v-x*((adj(x)*v)[0,0]/(adj(x)*x)[0,0].real)
    def sine(delta,v):
        lower=norm(v).lower()-delta
        assert lower>0
        return (delta/lower).upper()
    free_difference=max(sine(dw,orthogonal(s,w0)),sine(dz,orthogonal(t,z0)))
    old_difference=2*max(sine(dw,s-w0),sine(dz,t-z0))
    u0=free0*seed
    gram=adj(u0)*u0
    positive(gram)
    def edge(g,i,sign):
        return g[i,i].real+sign*sum(abs(g[i,j]) for j in range(g.ncols()) if j!=i)
    min_singular=min(edge(gram,i,-1).lower() for i in range(6)).sqrt()
    fullgram=adj(seed)*seed
    max_seed=max(edge(fullgram,i,1).upper() for i in range(6)).sqrt()
    error=free_difference*max_seed
    assert min_singular>error
    deviation=old_difference+2*error/(min_singular-error)
    mid_x=old0-2*span(u0)
    down=ranktwo['read'](before,after,weights,mid_x)
    positive(arb('.0129')*ranktwo['eye'](4)-block([v.upper() for v in down],after))
    maximum=(arb('.0129').sqrt()+deviation*arb('.293784').sqrt())**2
    assert maximum<arb('.01392'),maximum
    # 上型高字符线上的上型本征值为最高的一个。下型该字符来源全部同权。
    omega=acb(-arb(1)/2,arb(3).sqrt()/2)
    capvec=acb_mat([[0],[1],[omega],[omega**2]])*(1/arb(3).sqrt())
    assert (raw_up*capvec-capvec*eigen[1]).contains(acb_mat(4,1))
    # 高端额外投影取源下型大字符空间的六个本征态。其谱权均相同。
    native_high=[j for j,v in enumerate(weights) if arb('.1537')<v<arb('.1538')]
    assert len(native_high)==6
    other_character=acb_mat([[before[i,j] for j in native_high] for i in range(19)])
    assert (adj(other_character)*s).contains(acb_mat(6,1))
    assert (adj(other_character)*t).contains(acb_mat(6,1))
    assert (adj(other_character)*joint['align'](col(T,2),t)).contains(acb_mat(6,1))
    hi=ranktwo['read'](before,after,weights,old-2*other_character*adj(other_character))
    rawhi=block(hi,after)
    char=(adj(capvec)*rawhi*capvec)[0,0].real
    assert char>arb('.153')
    floor=(1-probs[1])*char
    assert floor>arb('.138') and maximum<arb('.01392')<arb('.01410')<floor
    # 下界适用于在固定六维核心上同作用的任意酉补全。只需平凡部门的四列。
    # P 的原始两列比从十九个站位列中挑选更稳定。
    seed_up=acb_mat([[acb(i,j) for i,j in zip(rr,ii,strict=True)]
                     for rr,ii in zip(ranktwo['REAL'],ranktwo['IMAG'],strict=True)])
    kp=E*seed_up*(arb(1)/10000)
    def four(p1,p2):
        return acb_mat([[kp[i,0],kp[i,1],p1[i,0],p2[i,0]] for i in range(19)])
    cols4=four(s,w)
    image4=four(w,s)
    for i in range(19):
        for j in range(2):
            image4[i,j]*=acb(0,1)
    native_trivial=[j for j in range(19) if sum(abs(v)**2 for v in (adj(E)*col(after,j)).entries())>arb('.5')]
    assert native_trivial==[0,5,6,9,11,15,18]
    gramv=adj(before)*cols4
    lower=[arb(0)]*19
    for j,shift in zip(native_trivial,('.00017','.0013','.0128','.0208','.00041','.0023','.082'),strict=True):
        gamma=arb(shift)
        rhs=adj(image4)*col(after,j)
        M=adj(gramv)*acb_mat([[gramv[i,k]/(weights[i]+gamma) for k in range(4)] for i in range(19)])
        # 以近似解挑选一个确切的试探向量。下界本身用实区间核验，不要求逆的误差估计。
        test=M.mid().solve(rhs.mid())
        test=test.mid()
        dual=2*(adj(test)*rhs)[0,0].real-(adj(test)*M*test)[0,0].real-gamma
        lower[j]=max(arb(0),dual.lower())
    trivial=acb_mat([[1,0],[0,1/arb(3).sqrt()],[0,1/arb(3).sqrt()],[0,1/arb(3).sqrt()]])
    bound_block=adj(trivial)*block(lower,after)*trivial
    level=arb('.0035')
    positive(bound_block-level*ranktwo['eye'](2))
    strange=level*(probs[1]+probs[3])
    assert strange>arb('.00064')
    print('全部自由酉补全的奇夸克下界',strange.lower().str(12),flush=True)
    print('低端点四阶上界、高端点每个上型中心三平面的底夸克下界',
          maximum.str(6),floor.lower().str(12),flush=True)


if __name__=='__main__':
    verify(Path(__file__).resolve().parent)
