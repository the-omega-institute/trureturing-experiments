"""实际十九站固定谱的正权振幅中心见证与指定平方质量见证的不可达性。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def adj(m):
    return m.conjugate().transpose()


def column(xs):
    return acb_mat([[x] for x in xs])


def norm(v):
    return (adj(v)*v)[0, 0].real.sqrt()


def coord(t):
    return [t[0, 0].real, t[1, 1].real, t[0, 1].real, t[0, 1].imag,
            t[2, 2].real, t[2, 3].imag]


def target_basis(up, down):
    vals, vec = (up+down).eig(right=True, algorithm="rump")
    seq = sorted(range(3), key=lambda j: vals[j].real.mid())
    assert all(vals[i].real < vals[j].real for i,j in zip(seq, seq[1:]))
    y = column([vec[i,seq[-1]] for i in range(3)])
    y = y * (y[0,0].conjugate()/(abs(y[0,0])*norm(y)))
    z = up*y - y*(adj(y)*up*y)[0,0].real
    z = z*(acb(0,1)/norm(z))
    x = column([(y[(i+1)%3,0]*z[(i+2)%3,0] -
                 y[(i+2)%3,0]*z[(i+1)%3,0]).conjugate() for i in range(3)])
    return acb_mat([[v[i,0] for v in (x,y,z)] for i in range(3)])


def native(root, number=2_306_888_100_085):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 320
    sites = shell["geo"]["sites"]
    rot = [sites.index((y,z,x)) for x,y,z in sites]
    assert (19+2*sum(i==j for i,j in enumerate(rot)))//3 == 7
    keep = [shell["soft"]["pole"], *shell["boundary_sites"]()]
    assert keep == [14,8,7,6]
    s,t,u = [1/arb(n).sqrt() for n in (3,2,6)]
    basis = acb_mat([[1,0,0,0], [0,s,t,u], [0,s,-t,u], [0,s,0,-2*u]])
    height = acb.zeta_zero(number).imag
    sectors = []
    for charge in (-4,2):
        h = shell["laplacian"](charge,height)
        assert all(h[i,j].overlaps(h[rot[i],rot[j]]) for i in range(19) for j in range(19))
        vals, vec = h.eig(right=True, algorithm="rump")
        seq = sorted(range(19), key=lambda j: vals[j].real.mid())
        assert all(vals[i].real < vals[j].real for i,j in zip(seq,seq[1:]))
        rows = []
        for j in seq:
            v = adj(basis)*column([vec[i,j] for i in keep])
            rows.append(v*adj(v)*(1/sum(abs(vec[i,j])**2 for i in range(19))))
        sectors.append(rows)
    return sectors


def bloch(block):
    trace = block[0,0].real+block[1,1].real
    return [2*block[0,1].real/trace, -2*block[0,1].imag/trace,
            (block[0,0].real-block[1,1].real)/trace]


def physical_frame(a,angle):
    sn,cs = angle.sin_cos()
    b = (1-a*a).sqrt()
    return acb_mat([[cs,0,-b*sn], [sn,0,b*cs], [0,1,0], [0,0,-a]])


def decompose(rows, m, a, angle, indices, *, return_weights=False):
    p,s,t = (m[i,i].real for i in range(3))
    k,w = m[0,2],m[1,2].imag
    b = (1-a*a).sqrt()
    gain = acb(0,arb("1.064")).exp()
    sn,cs = angle.sin_cos()
    orient = acb_mat([[cs,-sn,0,0],[sn,cs,0,0],[0,0,1,0],[0,0,0,1]])
    k = k*gain.conjugate()
    block = orient*acb_mat([[p,k/b,0,0],[k.conjugate()/b,(t-a*a*s)/(b*b),0,0],
                            [0,0,s,acb(0,-w/a)],[0,0,acb(0,w/a),s]])*adj(orient)
    picks = [i-1 for i in indices]
    helicities = [(v[2,3].imag.mid(),i) for i,v in enumerate(rows)]
    picks.extend((min(helicities)[1],max(helicities)[1]))
    assert len(set(picks)) == 6
    matrix = arb_mat([coord(rows[i]) for i in picks]).transpose()
    assert not matrix.det().contains(0)
    floor = arb("1e-9")
    identity = acb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    rhs = arb_mat([[v] for v in coord(block-identity*floor)])
    add = matrix.solve(rhs)
    weights = [floor for _ in rows]
    for i,wgt in zip(picks,add.entries(),strict=True):
        assert wgt > 0
        weights[i] += wgt
    assembled = sum((v*coef for v,coef in zip(rows,weights)),acb_mat(4,4))
    assert (assembled-block).contains(acb_mat(4,4))
    print("额外响应的一基谱序号及归一化增量",[(i+1,add[j,0].str(10)) for j,i in enumerate(picks)])
    return (assembled, weights) if return_weights else assembled


def dot(v,w):
    return sum(x*y for x,y in zip(v,w,strict=True))


def facets(rows):
    rays = [bloch(r) for r in rows if r[0,0].real+r[1,1].real>0]
    assert len(rays) == 7
    faces = []
    for a,b,c in combinations(rays,3):
        x,y = [[u-v for u,v in zip(p,a)] for p in (b,c)]
        n = [x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3] for i in range(3)]
        length = dot(n,n).sqrt()
        for sign in (-1,1):
            unit = [v*sign/length for v in n]
            off = max(dot(unit,r).upper() for r in rays)
            if (off-dot(unit,a)).upper()<arb("1e-20"):
                faces.append((unit,off))
    assert len(faces) == 10
    return faces


def exclude_pole_frame_amplitude(blocks, up_rows):
    # x 为极点时，原生平凡块的第一个对角元与相位无关。
    p,s,t = [blocks[0][j,j].real for j in range(3)]
    assert s > t
    upper = min(((abs(m[1,2].imag)/m[1,1].real)**2 +
                 m.det().real/(m[0,0].real*m[1,1].real**2)).sqrt().upper() for m in blocks)
    upper_fraction = p/(p+(t-upper**2*s)/(1-upper**2))
    fractions = [r[0,0].real/(r[0,0].real+r[1,1].real) for r in up_rows
                 if r[0,0].real+r[1,1].real>0]
    assert len(fractions) == 7
    floor = min(v.lower() for v in fractions)
    assert upper_fraction < arb(".02") < floor
    print("极点为 x 的一次幂见证所需极点份额上界与七个原生方向的下界",
          upper_fraction.str(11),floor.str(11))


def exclude_fixed_square(blocks, sectors):
    # 包含指定味标架的一切可行 a。正块的闭条件能同时处理半正定边缘。
    intervals = []
    for m in blocks:
        p,s = m[0,0].real,m[1,1].real
        lo = abs(m[1,2].imag)/s
        hi = (lo**2+m.det().real/(p*s*s)).sqrt()
        intervals.append((lo,hi))
    a = max(l.lower() for l,h in intervals).union(min(h.upper() for l,h in intervals))
    assert arb(".02145729")<a<arb(".02145745")
    b = (1-a*a).sqrt()
    vectors = []
    for m in blocks:
        p,s,t = [m[j,j].real for j in range(3)]
        k = m[0,2]/b
        vectors.append(bloch(acb_mat([[p,k],[k.conjugate(),(t-a*a*s)/(b*b)]])))
    faces = [facets(rows) for rows in sectors]
    # 对任意允许的 x 相位及平凡实基底转角，Bloch 向量依次绕 z、y 转动。
    stack = [(0,0,0)]
    tested,maxdepth = 0,0
    while stack:
        level,i,j = stack.pop()
        tested += 1
        maxdepth = max(maxdepth,level)
        unit = arb.pi()/(2**level)
        st,ct = ((2*i+1)*unit).sin_cos()
        sp,cp = ((2*j+1)*unit).sin_cos()
        # 每个半宽内移动的欧氏弦长不超过旋转角，正且迹一的 Bloch 半径至多一。
        pad = 2*unit
        separated = False
        for (x,y,z),domain in zip(vectors,faces,strict=True):
            u,v = cp*x-sp*y,sp*x+cp*y
            direction = (ct*u+st*z,v,-st*u+ct*z)
            if any(dot(n,direction)>limit+pad for n,limit in domain):
                separated = True
                break
        if not separated:
            assert level < 12
            stack.extend((level+1,2*i+di,2*j+dj) for di in (0,1) for dj in (0,1))
    print("平方质量同一重本征线的全相位全实转角排除盒数与深度",tested,maxdepth)


def verify(root):
    ckm = run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    spec = run_path(str(root/"ch39-three-gram-source-flavor.py"))["projections"]
    sectors = native(root)
    ctx.prec = 320
    v = acb_mat(ckm["standard_matrix"](*[arb(x).sin() for x in (".22704",".04275",".00377")],arb("1.139")))
    masses = [list(map(arb,row)) for row in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    a,angle = arb(".02145737"),arb("-.44525")
    frame = physical_frame(a,angle)
    for exponent in (1,2):
        diag = [acb_mat([[((m[i]/m[2])**exponent if i==j else 0) for j in range(3)]
                         for i in range(3)]) for m in masses]
        up,down = adj(v)*diag[0]*v,diag[1]
        q = target_basis(up,down)
        targets = [adj(q)*raw*q for raw in (up,down)]
        assert all(m[0,1].contains(0) and m[1,2].real.contains(0) for m in targets)
        if exponent == 2:
            exclude_fixed_square(targets,sectors)
            continue
        exclude_pole_frame_amplitude(targets,sectors[0])
        projs = []
        for rows,m,mass,picks in zip(sectors,targets,masses,((1,6,13,16),(1,6,10,16)),strict=True):
            f = decompose(rows,m,a,angle,picks)
            result = (adj(frame)*f*frame)*mass[-1]
            levels,pr = spec(result*result)
            assert all(l.overlaps(e*e) for l,e in zip(levels,mass,strict=True))
            projs.append(pr)
        pu,pd = projs
        assert all((pu[i]*pd[j]).trace().real.overlaps(abs(v[i,j])**2) for i in range(3) for j in range(3))
        j = (pu[0]*pd[0]*pu[1]*pd[1]).trace().imag
        assert j>0 and j.overlaps((v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag)
        print("实际十九站正振幅压缩的有符号 J",j.str(15))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
