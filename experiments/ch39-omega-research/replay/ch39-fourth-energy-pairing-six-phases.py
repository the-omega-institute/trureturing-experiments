"""N+3→N+4 全能序荷六配对，在首算术源定相上调整六个相位，精确实现 K、W、D。"""
from pathlib import Path
from runpy import run_path
from itertools import combinations

from flint import acb, acb_mat, arb, arb_mat, ctx


SIX = (17,15,0,3,18,9)
START = ("5.08216505","2.04530880","15.52360008","4.43705219","1.90924669","-2.55537555")


def adj(matrix):
    return matrix.conjugate().transpose()


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    second = lib("next-common-ckm-gravity")
    cx,rx,cy,ry,_ = second["verify"](root)
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    ordered = lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=n+1)[2][1]
    assert all(q is ordered[0] for q in ordered[:17])
    low,h,k = [ordered[j].sqrt() for j in (0,17,18)]
    native = lib("native-rh-ray-pvm")
    _,_,cross,_ = next(r for r in native["shell_weights"](root,n+4) if r[0]==6)
    c = [v.real for v in cross[0]]
    moments,heat,fields = lib("three-index-free-atoms")["broken_fields"](root,[n+4])
    field = fields[0]
    ctx.prec = 384
    weak = native["weak_series"](56)
    high = lib("full-neutral-free-cutoff")["weak_column"]
    mz,mp = arb("91.1876"),arb("1.22089e19")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    upper_weak = arb("-.074256577194663365545")
    near,far = arb("1.28"),arb(2)
    def rough(t):
        return c[15]*weak(t)+(1-c[15])*upper_weak
    assert rough(near)>0>rough(far)
    for _ in range(42):
        mid = (near+far)/2
        if rough(mid)>0:
            near = mid
        else:
            assert rough(mid)<0
            far = mid
    near -= arb("1e-11")
    far += arb("1e-11")
    time = near.union(far)
    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"]
    local = lib("next-common-ckm")["local"]
    def generator(center,radius,name,directions):
        box = [v+arb(0,radius.upper()) for v in center]
        base = seed(root,rank["eye"],name)
        return local(base,box,directions,rank["eye"](19),False)[0]
    x = generator(cx,rx,"ch39-full-ckm-gravity-cayley.tsv",lib("one-generator-gravity-full-ckm")["DIRECTIONS"])
    y = generator(cy,ry,"ch39-next-common-ckm-gravity-cayley.tsv",second["DIRECTIONS"])
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    source = spec(shell,6,n+1)[1]
    source = acb_mat([[source[i,17],source[i,18]] for i in range(19)])
    gram = native["gram_cholesky"](root)
    sites = (14,8,7,6)
    ends = [spec(shell,6,node)[1] for node in (n+3,n+4)]
    def frame(vectors,coeff):
        state = acb_mat([[coeff[sites.index(i)] if i in sites else 0] for i in range(19)])
        state *= 1/(adj(state)*state)[0,0].real.sqrt()
        ov = (adj(vectors)*state).entries()
        assert min(map(abs,ov))>arb(".00016")
        phases = [z/abs(z) for z in ov]
        return acb_mat([[vectors[i,j]*phases[j] for j in range(19)] for i in range(19)])
    s,t = [frame(vectors,gram[1]) for vectors in ends]
    v = adj(s)*t
    incoming = adj(s)*(y**3)*(x**3)*source
    def raw(point,derivative):
        phases = [acb(1)]*19
        for i,a in zip(SIX,point,strict=True):
            phases[i] = acb(0,a).exp()
        # 取掉最后一个目的侧单位相位，不影响两列的逐行模方。
        vd = acb_mat([[v[i,j]*phases[j] for j in range(19)] for i in range(19)])
        first = vd*incoming
        out = vd*first
        derivatives = []
        if derivative:
            vv = vd*vd
            for l in SIX:
                derivatives.append(acb_mat([[acb(0,1)*(vd[i,l]*first[l,j]+vv[i,l]*incoming[l,j])
                                              for j in (0,1)] for i in range(19)]))
        return out,derivatives

    def amplitudes(ov):
        return [low+(h-low)*abs(ov[i,0])**2+(k-low)*abs(ov[i,1])**2 for i in range(19)]

    base = amplitudes(raw([arb(0)]*6,False)[0])
    baseline_thermal = [moments(field,heat,z*z) for z in base]
    cols = (*zip(*baseline_thermal),[high(z*z,p0,p1) for z in base])
    readings = [sum(v*b for v,b in zip(c,col,strict=True)) for col in cols]
    assert min(base)>arb(".002") and readings[2]<arb("-.07")
    print("算术定相的第16层、牛顿、真空和弱差",*[z.str(12) for z in (base[15],*readings)],flush=True)
    # 单次换位只置换 V 的列及末壳的输出编号。先按统一截止下界排除整族弱差。
    residuals = [[int(j==k) for j in range(4)] for k in (1,2,3)]
    prev = (y**3)*(x**3)*source
    swaps = [()]+list(combinations(range(19),2))
    assert len(swaps)==172
    floors = []
    for coeff in gram[1:]+residuals:
        inside,outside = [frame(vectors,coeff) for vectors in ends]
        cross = adj(inside)*outside
        cin = adj(inside)*prev
        floor = arb(1)
        for swap in swaps:
            order = list(range(19))
            if swap:
                i,j = swap
                order[i],order[j] = order[j],order[i]
            permuted = acb_mat([[cross[i,j] for j in order] for i in range(19)])
            values = amplitudes(permuted*(permuted*cin))
            floor = min(floor,min(z.lower() for z in values))
        floors.append(floor)
    assert min(floors)>arb(".00042")
    # 高截止核在此段始终位于此常数加对称余量内。交叉权的总和严格为一。
    cutoff_floor = min(floors)**2
    assert sum(c).contains(1)
    weak_error = arb("4e-21")+arb("7.468e-36")/cutoff_floor
    weak_range = upper_weak+arb(0,(weak_error*sum(map(abs,c))).upper())
    assert weak_range<arb("-.07425")
    print("六种算术源在172种置换上的最小振幅下界及弱差统一包围",
          [z.str(11) for z in floors],weak_range.str(14),flush=True)

    def system(point,time):
        ov,dov = raw(point,True)
        f = amplitudes(ov)
        cutoff = mz*time.sqrt()/mp
        f[15] = cutoff
        eta = ((cutoff-low)/(h-low)).sqrt()
        residual = ov[15,0]-eta
        val = [residual.real,residual.imag,ov[15,1].real,ov[15,1].imag]
        jac = [[getattr(d[15,j],part) for d in dov] for j in (0,1) for part in ("real","imag")]
        err = [abs(c[15])*cutoff**2*sum(abs(a) for _,a,_ in field)/(2*arb.pi()),
               abs(c[15])*cutoff**4*sum(abs(g) for g,_,_ in field)/(8*arb.pi())]
        thermal = [arb(0,z.upper())-ref for z,ref in zip(err,(arb(1),arb("1e-122")),strict=True)]
        slopes = [[],[]]
        for i,z in enumerate(f):
            if i==15:
                for row in slopes:
                    row.append(arb(0))
                continue
            assert z>arb("1e-6")
            for r,a in enumerate(moments(field,heat,z*z)):
                thermal[r] += c[i]*a
            ex = [(-u/z**2).exp() for _,_,u in field]
            slopes[0].append(c[i]*z*sum(a*e for (_,a,_),e in zip(field,ex,strict=True))/arb.pi())
            slopes[1].append(-c[i]*z**3*sum(g*e for (g,_,_),e in zip(field,ex,strict=True))/(2*arb.pi()))
        df = arb_mat([[2*((h-low)*(ov[i,0].conjugate()*d[i,0]).real+
                          (k-low)*(ov[i,1].conjugate()*d[i,1]).real) for d in dov] for i in range(19)])
        jac.extend((arb_mat(slopes)*df).tolist())
        val.extend(thermal)
        return arb_mat([[z] for z in val]),arb_mat(jac),f

    center = list(map(arb,START))
    for _ in range(4):
        val,jac,_ = system(center,time.mid())
        step = jac.mid().solve(val.mid())
        center = [(v-w).mid() for v,w in zip(center,step.entries(),strict=True)]
    val,jac,_ = system(center,time)
    pre = jac.mid().inv()
    assert not pre.det().contains(0)
    radius = arb("1e-10")
    box = [v+arb(0,radius.upper()) for v in center]
    _,wide,amps = system(box,time)
    defect = arb_mat([[int(i==j) for j in range(6)] for i in range(6)])-pre*wide
    lip = max(sum(abs(z).upper() for z in row) for row in defect.tolist())
    offset = max(abs(z).upper() for z in (pre*val).entries())/radius
    print("六相位中心、盒半宽、位移占比和导数界",
          [z.str(14) for z in center],radius.str(5),offset.str(12),lip.str(12),flush=True)
    assert lip+offset<1
    other = sum(c[i]*high(z*z,p0,p1) for i,z in enumerate(amps) if i!=15)
    ends = [c[15]*weak(t)+other for t in (near,far)]
    assert ends[0]>0>ends[1]
    print("第16层交叉权、截止 GeV、两个弱差符号",
          *[z.str(13) for z in (c[15],mz*time.sqrt(),*ends)],flush=True)
    print("第18、19层及其余十八层的振幅下界",
          amps[17].str(12),amps[18].str(12),min(z.lower() for i,z in enumerate(amps) if i!=15).str(12),flush=True)
    return center,radius,time


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
