"""右侧相邻壳同一十九维物理酉识别的三荷最低层 Gram 判据与见证。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


CENTER = 2_306_888_100_085


def adj(x):
    return x.conjugate().transpose()


def cols(a, indices):
    return acb_mat([[a[i,j] for j in indices] for i in range(a.nrows())])


def join(*matrices):
    return acb_mat([[x[i,j] for x in matrices for j in range(x.ncols())]
                    for i in range(matrices[0].nrows())])


def eye(n):
    return acb_mat([[int(i == j) for j in range(n)] for i in range(n)])


def length(x):
    return (adj(x)*x)[0,0].real.sqrt()


def ortho_columns(columns, *, start=(), take=None):
    """投影到已选方向的正交补，按剩余范数的下界选择补基。"""
    chosen = list(start)
    extra = []
    count = len(columns) if take is None else take
    for _ in range(count):
        components = []
        for v in columns:
            projection = sum((a*(adj(a)*v)[0,0] for a in chosen), acb_mat(v.nrows(),1))
            components.append(v-projection)
        i = max(range(len(components)), key=lambda j: length(components[j]).lower())
        v = components[i]
        n = length(v)
        assert n > arb('1e-8')
        v = v*(1/n)
        chosen.append(v)
        extra.append(v)
        columns = columns[:i] + columns[i+1:]
    return extra


def orthogonal_pair(h, low):
    """源夸克三维最低空间中取正交于荷六第十八层的二维平面。"""
    r = adj(h)*low
    assert abs(r[0,2]) > arb('.1')
    ns = acb_mat([[1,0],[0,1],[-r[0,0]/r[0,2],-r[0,1]/r[0,2]]])
    candidates = low*ns
    return join(*ortho_columns([cols(candidates,[j]) for j in range(2)]))


def eigensystem(h):
    values,vectors = h.eig(right=True,algorithm='rump')
    order = sorted(range(len(values)),key=lambda i: values[i].real.mid())
    assert all(values[j].imag.contains(0) for j in order)
    vals = [values[j].real for j in order]
    assert all(a<b for a,b in zip(vals,vals[1:]))
    vecs = [cols(vectors,[j]) for j in order]
    return vals,[v*(1/length(v)) for v in vecs]


def spinor(u):
    a,b = u[0,0],u[1,0]
    return acb_mat([[a,-b.conjugate()],[b,a.conjugate()]])*(1/length(u))


def rotation_from_grams(t, y):
    """按共同列序作相同的正交化，再以站点坐标分别补齐两侧。"""
    def full(a):
        basis = []
        for j in range(a.ncols()):
            basis += ortho_columns([cols(a,[j])], start=basis, take=1)
        units = eye(19)
        basis += ortho_columns([cols(units,[j]) for j in range(19)], start=basis, take=19-len(basis))
        return join(*basis)
    u = full(t)*adj(full(y))
    error = max(abs(z).upper() for z in (adj(u)*t-y).entries())
    unit_error = max(abs(z).upper() for z in (adj(u)*u-eye(19)).entries())
    assert max(error,unit_error) < arb('1e-16'), (error,unit_error)
    return u


def candidate(source, target, q, pair, single):
    other = 2 if q == -4 else -4
    tq = cols(target[q],pair)
    to = cols(target[other],[single])
    g = cols(target[6],[0,15])
    l6 = cols(source[6], range(17))
    h = cols(source[6],[17])
    a,b = [orthogonal_pair(h,cols(source[charge],[0,15,17])) for charge in (q,other)]
    cross = adj(a)*b
    vals,vecs = eigensystem(adj(cross)*cross)
    c = adj(tq)*to
    eta = (length(c)**2-vals[0])/(vals[1]-vals[0])
    assert 0 < eta < 1
    z = b*(vecs[0]*(1-eta).sqrt()+vecs[1]*eta.sqrt())
    yq = a*spinor(adj(a)*z)*adj(spinor(c))
    Y,T = join(yq,z),join(tq,to)
    assert (adj(Y)*Y-adj(T)*T).contains(acb_mat(3,3))
    py = l6*(adj(l6)*Y)
    cross_g = adj(T)*g
    initial = py*(adj(py)*py).solve(cross_g)
    remainder = eye(2)-adj(initial)*initial
    eigen,_ = eigensystem(remainder)
    assert eigen[0] > arb('.0002')
    # 增加落在 span(PY) 正交补中的两列，补足荷六目标所需 Gram。
    pbase = ortho_columns([cols(py,[j]) for j in range(3)])
    axes = ortho_columns([cols(l6,[j]) for j in range(17)],start=pbase,take=2)
    d = remainder[0,0].real.sqrt()
    off = remainder[0,1]/d
    cholesky = acb_mat([[d,off],[0,(remainder[1,1].real-abs(off)**2).sqrt()]])
    yg = initial+join(*axes)*cholesky
    t,y = join(T,g),join(Y,yg)
    assert (adj(t)*t-adj(y)*y).contains(acb_mat(5,5))
    u = rotation_from_grams(t,y)
    print('双层荷与层、另一荷单层、交叉范数及其允许区间、荷六剩余 Gram 最小本征值',
          q,tuple(i+1 for i in pair),single+1,length(c).str(11),
          *[v.sqrt().str(11) for v in vals],eigen[0].str(11),flush=True)
    return u


def verify(root):
    shells = run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    pvm = run_path(str(root/'ch39-adjacent-pvm-transport.py'))
    three = run_path(str(root/'ch39-three-shell-two-space-readout.py'))
    train,_ = three['verify'](root,return_data=True)
    ctx.prec = 768
    source,target = {},{}
    for off,out in ((1,source),(2,target)):
        for q in (-4,2,6):
            levels,out[q] = pvm['spectrum'](shells,q,CENTER+off)
            if off == 1:
                row = pvm['training_values'](levels,train)
                low = range(17) if q == 6 else (0,15,17)
                base = row[0]
                assert all(row[i].overlaps(base) if i in low else row[i] > 2*base for i in range(19))
                if q != 6:
                    assert base.contains(arb('1e-9'))
    ctx.prec = 768
    # 全取一侧三层时，另一侧单向量投影范数必落在两源平面的奇异值区间。
    s = adj(cols(source[-4],[0,15,17]))*cols(source[2],[0,15,17])
    src,_ = eigensystem(adj(s)*s)
    assert src[0].sqrt() > arb('.9782')
    for q,triple,single in ((-4,(0,15,17),15),(2,(0,15,18),15)):
        other = 2 if q == -4 else -4
        n = length(adj(cols(target[q],triple))*cols(target[other],[single]))
        assert n < arb('.96')
        e = arb('.0085') if q == -4 else arb('.0105')
        assert src[0].sqrt()*(1-e**2).sqrt() - (1+arb(3).sqrt())*e > n
        print('排除的三层荷、另一荷单层、必要源下界、目标投影范数',
              q,single+1,src[0].sqrt().str(12),n.str(12),'至少一层的泄漏概率大于', (e**2).str(8),flush=True)
    return {q:candidate(source,target,q,(0,15),0) for q in (-4,2)}


if __name__=='__main__':
    verify(Path(__file__).resolve().parent)
