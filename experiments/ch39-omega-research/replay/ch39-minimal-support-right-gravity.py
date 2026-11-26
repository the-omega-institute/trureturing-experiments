"""右侧单步荷六三项精确匹配的最小非平凡站位支撑。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, arb_mat, ctx


CENTER = 2_306_888_100_085

def eye(n):
    return acb_mat([[int(i == j) for j in range(n)] for i in range(n)])

def adj(M):
    return M.conjugate().transpose()

def column_norm(v):
    return sum(abs(z)**2 for z in v.entries()).sqrt()

def basis3():
    seeds = [
        [(' .29400263','.35913994'),('-.23345513','-.52693317'),('-.07035766','-.66893265')],
        [('.74840468','.15258358'),('.58048402','-.04496589'),('-.11134556','.25539620')],
        [('.44136486','.07991686'),('-.44906009','.35665325'),('.68432311','.04069071')],
    ]
    cols = []
    for col in zip(*seeds, strict=True):
        v = acb_mat([[acb(x.strip(),y)] for x,y in col])
        for u in cols:
            v -= u*(adj(u)*v)[0,0]
        cols.append(v*(1/column_norm(v)))
    return acb_mat([[col[j,0] for col in cols] for j in range(3)])

def charts():
    # Cayley 坐标的厄米生成元先列对角元，再按字典序列每对指标的实、虚元。
    def generators(n, edges):
        out = []
        for j in range(n):
            H = acb_mat(n,n); H[j,j] = 1; out.append(H)
        for i,j in edges:
            for z in (acb(1), acb(0,1)):
                H = acb_mat(n,n); H[i,j] = z; H[j,i] = z.conjugate(); out.append(H)
        return out
    g = generators(3,combinations(range(3),2))
    yield [14,6,7], basis3(), [g[k] for k in (0,2,3,6,7,8)]
    theta, a,b,c, *extras = map(arb, ('-.61467505','.24569355','-1.75750955','-.61541768',
                                         '-.59978879','-2.97722605','-1.50854174','2.89052335'))
    def cis(x): return acb(0,x).exp()
    s,h = theta.sin_cos()
    U = acb_mat(6,6)
    U[0,0],U[0,1] = cis(a)*h,cis(b)*s
    U[1,0],U[1,1] = -cis(c-b)*s,cis(c-a)*h
    for i,x in enumerate(extras,2): U[i,i] = cis(x)
    g = generators(6,[(0,1)])
    yield [14,17,8,16,4,13],U,[g[k] for k in (0,1,2,3,4,7)]


def fixed_two(source, target, high, mz_ratio):
    overlaps = (adj(target)*source)
    gaps = []
    for site in range(19):
        if site == 14: continue
        sites = [14, site]
        A = acb_mat([[source[j,k] for k in (17,18)] for j in sites])
        y = acb_mat([[target[j,15].conjugate() for j in sites]])
        tail = acb_mat([[overlaps[15,k] for k in (17,18)]])-y*A
        offset = tail*A.inv()
        # σ_min(A)≥|det A|/||A||_F。此处只需严格正的统一下界。
        gap = abs(column_norm(offset)-column_norm(y))*abs(A.det())/column_norm(A)
        gaps.append(gap.lower())
    gap = min(gaps)
    assert len(gaps) == 18 and high*gap**2 > 16*mz_ratio
    print('外站固定的十八种两站块的重叠向量模及第16层振幅下界', gap.str(11),(high*gap**2).lower().str(11),flush=True)


def witness(sites, anchor, generators, source, target, levels, weights, field, moments, heat, low_weak, high_weak):
    n = len(sites); I = eye(n); im = acb(0,1)
    assert (adj(anchor)*anchor-I).contains(acb_mat(n,n))
    A = acb_mat([[source[j,k] for k in (17,18)] for j in sites])
    B = acb_mat([[target[j,k].conjugate() for j in sites] for k in range(19)])
    full = adj(target)*source
    rest = acb_mat([[full[j,k] for k in (17,18)] for j in range(19)]) - B*A
    low, m, h = levels
    pi = arb.pi(); p0=(arb('91.1876')/arb('1.22089e19'))**2; p1=(arb('1.176397e14')/arb('1.22089e19'))**2
    t0, t1 = arb('1.30'), arb('1.50')
    reference = arb('-.074256577194663365545')
    def rough(t): return weights[15]*low_weak(t)+(1-weights[15])*reference
    assert rough(t0)>0>rough(t1)
    for _ in range(20):
        z=(t0+t1)/2
        if rough(z)>0:t0=z
        else:
            assert rough(z)<0
            t1=z
    def system(point,t,*,derivative_only=False):
        special = (p0*t).sqrt()
        eta = ((special-low)/(m-low)).sqrt()
        H = sum((G*x for G,x in zip(generators,point,strict=True)),acb_mat(n,n))
        inv=(I-im*H).inv()
        U=(2*inv-I)*anchor
        M=rest+B*U*A
        dM=[B*inv*G*inv*anchor*A*(2*im) for G in generators]
        vals=[M[15,0].real-eta,M[15,1].real,M[15,0].imag,M[15,1].imag,-arb(1),-arb('1e-122')]
        jac=arb_mat(6,6)
        for k, ds in enumerate(dM):
            for j,val in enumerate((ds[15,0].real,ds[15,1].real,ds[15,0].imag,ds[15,1].imag)):jac[j,k]=val
        fixed_weak=arb(0);others=[]
        for j,cj in enumerate(weights):
            F=special if j==15 else low+(m-low)*abs(M[j,0])**2+(h-low)*abs(M[j,1])**2
            assert F>0
            x=F**2
            if not derivative_only:
                K,W=moments(field,heat,x);vals[4]+=cj*K;vals[5]+=cj*W
            if j==15:continue
            fixed_weak+=cj*high_weak(x,p0,p1);others.append(F)
            decay=[(-z/x).exp() for _,_,z in field]
            dK=cj*F*sum(a*e for (_,a,_),e in zip(field,decay,strict=True))/pi
            dW=-cj*F**3*sum(g*e for (g,_,_),e in zip(field,decay,strict=True))/(2*pi)
            for k,ds in enumerate(dM):
                dF=2*((m-low)*(M[j,0].conjugate()*ds[j,0]).real+(h-low)*(M[j,1].conjugate()*ds[j,1]).real)
                jac[4,k]+=dK*dF;jac[5,k]+=dW*dF
        return arb_mat([[v] for v in vals]),jac,fixed_weak,min(v.lower() for v in others)
    point=[arb(0)]*6
    for _ in range(5):
        F,J,_,_=system(point,(t0+t1)/2)
        point=[(v-d).mid() for v,d in zip(point,J.solve(F).entries(),strict=True)]
    radius=arb('1e-7');box=[x+arb(0,radius.upper()) for x in point];t=t0.union(t1)
    F,J,_,_=system(point,t);pre=J.inv();assert not pre.det().contains(0)
    _,wide,other,minimum=system(box,t,derivative_only=True)
    defect=arb_mat([[int(i==j) for j in range(6)] for i in range(6)])-pre*wide
    frac=max((abs(v)/radius+sum(abs(x) for x in row)).upper() for v,row in zip((pre*F).entries(),defect.tolist(),strict=True))
    ends=[weights[15]*low_weak(z)+other for z in (t0,t1)]
    print('局部站、六面误差占半宽、其余振幅下界、弱差两端', sites,frac.str(9),minimum.str(9),*[v.str(9) for v in ends],flush=True)
    assert frac<1 and ends[0]>0>ends[1]


def verify(root):
    sh=run_path(str(root/'ch39-three-ray-chiral-boundary.py'));trans=run_path(str(root/'ch39-adjacent-pvm-transport.py'))
    pv=run_path(str(root/'ch39-native-rh-ray-pvm.py'));gravity=run_path(str(root/'ch39-three-index-free-atoms.py'))
    high_weak=run_path(str(root/'ch39-full-neutral-free-cutoff.py'))['weak_column']
    ordered=run_path(str(root/'ch39-adjacent-ordered-charge-six-cutoffs.py'))['verify'](root)[2]
    _,pos=ordered; assert all(x is pos[0] for x in pos[:17]); levels=[pos[i].sqrt() for i in (0,17,18)]
    n=CENTER+2;moments,heat,fields=gravity['broken_fields'](root,[n]);ctx.prec=320
    _,src=trans['spectrum'](sh,6,n-1);_,dest=trans['spectrum'](sh,6,n)
    _,_,cross,_=next(t for t in pv['shell_weights'](root,n) if t[0]==6);ctx.prec=320
    weights=[z.real for z in cross[0]];assert sum(weights).contains(1)
    mz=arb('91.1876')/arb('1.22089e19');fixed_two(src,dest,levels[1],mz)
    weak=pv['weak_series']();ranges=run_path(str(root/'ch39-outer-shared-distance-kernel.py'))['weak_bands']
    whole,large=ranges(weak,high_weak,mz**2,(arb('1.176397e14')/arb('1.22089e19'))**2)
    assert sum(c*(large if i==15 else whole) for i,c in enumerate(weights))<arb('-.01656')
    for block,base,G in charts():witness(block,base,G,src,dest,levels,weights,fields[0],moments,heat,weak,high_weak)


if __name__=='__main__':verify(Path(__file__).resolve().parent)
