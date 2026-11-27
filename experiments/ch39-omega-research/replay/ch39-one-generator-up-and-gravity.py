"""右外壳以三种本征相位的同一发生元实现上型三质量与有轻场破缺的荷六三读数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def adj(a):
    return a.conjugate().transpose()


def inner(a,b):
    return (adj(a)*b)[0,0]


def col(a,j):
    return acb_mat([[a[i,j]] for i in range(a.nrows())])


def side_by_side(vectors):
    return acb_mat([[v[i,0] for v in vectors] for i in range(vectors[0].nrows())])


def align(v,s):
    q = inner(v,s)
    assert abs(q)>0
    return v*(q/abs(q))


def mode_parameters(P,E,S,T):
    s = col(S,17)
    R = E*adj(E)-P
    target = [align(col(T,k),s) for k in (15,17)]
    imag = acb(0,1)
    M = R*side_by_side([target[0],imag*target[0],target[1],imag*target[1],imag*s])
    gram = adj(M)*M
    gram_real = arb_mat([[v.real for v in row] for row in gram.tolist()])
    assert gram_real.det()>0
    dual = M*acb_mat(gram_real.inv())
    fixed = -imag*P*s
    vals = [-inner(z,fixed) for z in target]
    offset = dual*acb_mat([[x] for z in vals for x in (z.real,z.imag)]+[[0]])
    seed = R*align(col(T,0),s)
    tangent = seed-dual*acb_mat([[v.real] for v in (adj(M)*seed).entries()])
    tangent *= 1/inner(tangent,tangent).real.sqrt()
    return (fixed,offset,col(dual,0),col(dual,2),tangent,inner(R*s,R*s).real)


def verify(root, *, return_parameters=False):
    base = run_path(str(root/'ch39-rank-two-generator-up-quark.py'))
    orbits = run_path(str(root/'ch39-up-quark-arbitrary-site-phase-floor.py'))['ORBITS']
    shells = run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    transport = run_path(str(root/'ch39-adjacent-pvm-transport.py'))
    grav = run_path(str(root/'ch39-three-index-free-atoms.py'))
    pivots = run_path(str(root/'ch39-native-rh-ray-pvm.py'))
    high_weak = run_path(str(root/'ch39-full-neutral-free-cutoff.py'))['weak_column']
    n = base['N']+2
    _,positions = run_path(str(root/'ch39-adjacent-ordered-charge-six-cutoffs.py'))['verify'](root)[2]
    assert all(x is positions[0] for x in positions[:17])
    low,h,k = [positions[j].sqrt() for j in (0,17,18)]
    moments,heat,fields = grav['broken_fields'](root,[n])
    field = fields[0]
    _,_,cross,_ = next(x for x in pivots['shell_weights'](root,n) if x[0]==6)
    c = [z.real for z in cross[0]]
    ctx.prec = 384
    weak = pivots['weak_series']()
    S = transport['spectrum'](shells,6,n-1)[1]
    T = transport['spectrum'](shells,6,n)[1]
    E = acb_mat([[1/arb(len(o)).sqrt() if i in o else 0 for o in orbits] for i in range(19)])
    P = base['projector'](orbits)
    fixed,offset,small,big,d,rho = mode_parameters(P,E,S,T)
    trivial = [i for i in range(19) if inner(adj(E)*col(T,i),adj(E)*col(T,i)).real > arb('.5')]
    assert trivial==[0,3,6,9,12,15,17]
    # 简单谱和三循环对称性保证每条本征线具有唯一字符。
    assert inner(adj(E)*col(S,17),adj(E)*col(S,17)).real > arb('.99')
    assert inner(adj(E)*col(S,18),adj(E)*col(S,18)).real < arb('.01')
    s19 = col(S,18)
    char3,char19 = [align(col(T,j),s19) for j in (2,18)]
    assert all(abs(inner(s19,v))>0 for v in (char3,char19))
    high = [j for j in trivial if j!=15]+[2,18]
    quiet = sum(v for j,v in enumerate(c) if j not in high and j!=15)
    mz,mp = arb('91.1876'),arb('1.22089e19')
    p0,p1 = (mz/mp)**2,(arb('1.176397e14')/mp)**2

    def vectors(a,b,gev):
        e = ((gev/mp-low)/(h-low)).sqrt()
        linear = offset+small*e+big*a
        gap = rho-inner(linear,linear).real
        assert gap>0
        y = fixed+linear+d*gap.sqrt()
        yp = big-d*(inner(linear,big).real/gap.sqrt())
        z = char19*b+char3*(1-b*b).sqrt()
        return y,yp,z

    def amplitudes(a,b,gev):
        y,yp,_ = vectors(a,b,gev)
        f = [low]*19
        da,db = [arb(0)]*19,[arb(0)]*19
        f[15] = gev/mp
        for j in trivial:
            if j!=15:
                v = col(T,j)
                p = inner(v,y)
                f[j] = low+(h-low)*abs(p)**2
                da[j] = 2*(h-low)*(p.conjugate()*inner(v,yp)).real
        f[2] = low+(k-low)*(1-b*b)
        f[18] = low+(k-low)*b*b
        db[18] = 2*(k-low)*b
        db[2] = -db[18]
        return f,da,db

    def thermal(pair,gev):
        f,da,db = amplitudes(*pair,gev)
        # 所有低层质量用积分表示的绝对界外包，避开极低截止的数值相减。
        threshold = max(low.upper(),(gev/mp).upper())
        low_weight = abs(quiet)+abs(c[15])
        err = [low_weight*threshold**2*sum(abs(v) for _,v,_ in field)/(2*arb.pi()),
               low_weight*threshold**4*sum(abs(g) for g,_,_ in field)/(8*arb.pi())]
        val = [arb(0,e.upper())-goal for e,goal in zip(err,(1,arb('1e-122')),strict=True)]
        jac = arb_mat(2,2)
        for j in high:
            amp = f[j]
            assert amp>arb('.0001')
            x = amp**2
            cols = moments(field,heat,x)
            decay = [(-z/x).exp() for _,_,z in field]
            slopes = [amp*sum(v*e for (_,v,_),e in zip(field,decay,strict=True))/arb.pi(),
                      -amp**3*sum(g*e for (g,_,_),e in zip(field,decay,strict=True))/(2*arb.pi())]
            for i in range(2):
                val[i] += c[j]*cols[i]
                for l,dx in enumerate((da[j],db[j])):
                    jac[i,l] += c[j]*slopes[i]*dx
        return arb_mat([[v] for v in val]),jac

    l,r = arb(110),arb(115)
    pair=[arb('.859'),arb('.919')]
    for _ in range(5):
        val,J = thermal(pair,(l+r)/2)
        pair = [(x-y).mid() for x,y in zip(pair,J.solve(val).entries(),strict=True)]
    radius=arb('3e-6')
    box=[a+arb(0,radius.upper()) for a in pair]
    ref,_=thermal(pair,l.union(r))
    pre=thermal(pair,(l+r)/2)[1].inv()
    assert not pre.det().contains(0)
    mat=arb_mat([[1,0],[0,1]])-pre*thermal(box,l.union(r))[1]
    ratio=max((abs(v)/radius+sum(abs(x) for x in row)).upper()
              for v,row in zip((pre*ref).entries(),mat.tolist(),strict=True))
    assert ratio<1,ratio
    f,_,_=amplitudes(*box,l.union(r))
    other=quiet*weak(low**2/p0)+sum(c[j]*high_weak(f[j]**2,p0,p1) for j in high)
    def face(q):
        return other+c[15]*weak((q/mz)**2)
    assert face(l)>0>face(r)
    for _ in range(26):
        mid=(l+r)/2
        if face(mid)>0: l=mid
        elif face(mid)<0: r=mid
        else: break
    assert face(l)>0>face(r)
    q16=l.union(r)
    y,_,z=vectors(*box,q16)
    # 对两个正交字符分别构造把源高向量送到所求像的 Householder 反射。
    I=base['eye'](19)
    s=col(S,17)
    a=(I-P)*s-(y-fixed)
    b=s19-z
    na,nb=inner(a,a).real,inner(b,b).real
    assert na>0 and nb>0
    Q=a*adj(a)*(1/na)+b*adj(b)*(1/nb)
    assert (P*Q).contains(acb_mat(19,19)) and (adj(a)*b)[0,0].contains(0)
    X=I+(acb(0,1)-1)*P-2*Q
    assert (adj(X)*X-I).contains(acb_mat(19,19))
    assert (X*X-(I-2*P)).contains(acb_mat(19,19))
    image=X*X*X*side_by_side([s,s19])
    assert (image-side_by_side([y,z])).contains(acb_mat(19,2))
    print('两热参数、统一压缩比、弱差异号端点、第16层GeV',
          [x.str(12) for x in box],ratio.str(9),face(l).lower().str(8),face(r).upper().str(8),q16.str(13),flush=True)
    print('第18层与第19层振幅',f[17].str(9),f[18].str(9),flush=True)
    if return_parameters:
        return X, tuple(box), q16, (low,h,k)
    return X


if __name__=='__main__':
    verify(Path(__file__).resolve().parent)
