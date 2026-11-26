"""左侧荷六三项匹配的最小几何支撑，含全部 816 个四站禁区。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


CENTER = 2_306_888_100_085
ACTIVE = (6, 7, 17, 18)  # 零基谱列，对应一基第 7、8、18、19 层。
# 各项给出站位、Cayley 生成元序号、混合块的 Gram–Schmidt 种子、其余站的初相。
SEEDS = [([14, 11, 12, 18, 6], [5, 9, 10, 11, 12, 8, 7, 6, 18, 21],
  [[('0.5064229597', '-0.0924897700'), ('0.0158751958', '0.2361084684'), ('0.5118169489', '0.2760208038'),
    ('-0.1762128978', '-0.4304957197'), ('0.3522856532', '-0.0188518664')],
   [('0.0523016460', '-0.6177457332'), ('-0.2530060075', '-0.5508516062'), ('-0.2366086052', '0.4089355251'),
    ('0.0010594086', '-0.1122424140'), ('0.0241248985', '-0.1086822287')],
   [('0.0364047509', '0.3764029876'), ('-0.2903367529', '-0.2942400969'), ('-0.0584950212', '-0.3041453118'),
    ('-0.1528127041', '-0.5781078126'), ('-0.2256372396', '-0.4262912769')],
   [('0.3662193831', '0.0414136771'), ('-0.4130093915', '0.1660511140'), ('0.0125922943', '-0.0367412584'),
    ('0.3385894375', '0.4702045092'), ('0.1161410775', '-0.5615035328')],
   [('-0.1271129313', '-0.2358860613'), ('-0.2851490646', '0.3552729286'), ('-0.4512987189', '-0.3708624721'),
    ('0.0417061047', '-0.2756594198'), ('0.5329127982', '0.1331846603')]],
  []),
 ([14, 17, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 18], [20, 0, 10, 8, 19, 14, 13, 9, 15, 4],
  [[('-0.6198800674', '-0.2647040381'), ('0.3012194194', '-0.6744978396')],
   [('-0.7087534234', '0.2082043685'), ('0.1793300730', '0.6497386019')]],
  ['-2.9032848067', '-1.0932243661', '0.6882931364', '-1.4996099758', '-1.2918131746', '-1.4304245960', '-2.3513382270',
   '-2.8607789893', '-2.4246924854', '-1.3692085185', '-2.0784180271', '-2.6852959428', '0.7619104098', '-0.9870888383',
   '-0.4356330810', '-1.4431576374', '-1.0825191305'])]


def adj(matrix):
    return matrix.conjugate().transpose()


def identity(n):
    return acb_mat([[int(i == j) for j in range(n)] for i in range(n)])


def norm(matrix):
    return sum(abs(z)**2 for z in matrix.entries()).sqrt()


def chart(seed):
    sites, selected, raw, phases = seed
    n, m = len(sites), len(raw)
    cols = []
    for column in zip(*raw, strict=True):
        vec = acb_mat([[acb(x, y)] for x, y in column])
        for u in cols:
            vec -= u*(adj(u)*vec)[0,0]
        cols.append(vec*(1/norm(vec)))
    anchor = acb_mat(n,n)
    for i in range(m):
        for j in range(m):
            anchor[i,j] = cols[j][i,0]
    for j, phase in enumerate(phases, m):
        anchor[j,j] = acb(0, phase).exp()
    generators = []
    for i in range(n):
        h = acb_mat(n,n); h[i,i] = 1; generators.append(h)
    for i, j in combinations(range(m), 2):
        for z in (acb(1), acb(0,1)):
            h = acb_mat(n,n); h[i,j] = z; h[j,i] = z.conjugate(); generators.append(h)
    assert (adj(anchor)*anchor-identity(n)).contains(acb_mat(n,n))
    return sites, anchor, [generators[i] for i in selected]


def four_site_bound(source, target, mid, floor):
    full = adj(target)*source
    overlaps = acb_mat([[full[15,j] for j in ACTIVE]])
    lower = []
    for other in combinations([i for i in range(19) if i != 14], 3):
        row = (14, *other)
        A = acb_mat([[source[i,j] for j in ACTIVE] for i in row])
        y = acb_mat([[target[i,15].conjugate() for i in row]])
        assert not A.det().contains(0)
        inv = A.inv()
        b = overlaps - y*A
        lower.append((abs(norm(b*inv)-norm(y))/norm(inv)).lower())
    minimum = min(lower)
    assert len(lower) == 816 and mid*minimum**2 > 16*floor
    print('四站块个数、高层重叠行范数及第16层振幅的统一下界', len(lower),
          minimum.str(12), (mid*minimum**2).lower().str(12), flush=True)


def witness(seed, source, target, amplitudes, weights, field, moments, heat, weak, high_weak):
    sites, anchor, gs = chart(seed)
    n, ng = len(sites), len(gs)
    assert ng == 10
    I, im = identity(n), acb(0,1)
    A = acb_mat([[source[j,k] for k in ACTIVE] for j in sites])
    B = acb_mat([[target[j,k].conjugate() for j in sites] for k in range(19)])
    full = adj(target)*source
    rest = acb_mat([[full[j,k] for k in ACTIVE] for j in range(19)]) - B*A
    low, middle, high = amplitudes
    delta = [middle-low]*3 + [high-low]
    pi = arb.pi(); p0 = (arb('91.1876')/arb('1.22089e19'))**2
    p1 = (arb('1.176397e14')/arb('1.22089e19'))**2
    t0, t1 = arb('1.01'), arb('1.30')
    ref = arb('-.074256577194663365545')
    def rough(t):
        return weights[15]*weak(t)+(1-weights[15])*ref
    assert rough(t0) > 0 > rough(t1)
    for _ in range(20):
        z=(t0+t1)/2
        if rough(z)>0:t0=z
        else:
            assert rough(z)<0
            t1=z

    def system(point, t, *, derivative_only=False):
        special = (p0*t).sqrt()
        eta = ((special-low)/delta[0]).sqrt()
        H = sum((G*x for G,x in zip(gs,point,strict=True)),acb_mat(n,n))
        inv = (I-im*H).inv()
        M = rest+B*(2*inv-I)*anchor*A
        left = B*inv
        right = inv*anchor*A
        dM = [left*G*right*(2*im) for G in gs]
        vals = [M[15,k].real for k in range(4)]+[M[15,k].imag for k in range(4)]+[-arb(1),-arb('1e-122')]
        vals[0] -= eta
        J = arb_mat(ng,ng)
        for i, derivative in enumerate(dM):
            for k in range(4):
                J[k,i],J[k+4,i] = derivative[15,k].real,derivative[15,k].imag
        others = []; fixed_weak = arb(0)
        for j, c in enumerate(weights):
            F = special if j == 15 else low+sum(d*abs(M[j,k])**2 for k,d in enumerate(delta))
            assert F>0
            x = F**2
            if not derivative_only:
                kcol,wcol = moments(field,heat,x)
                vals[8] += c*kcol; vals[9] += c*wcol
            if j==15:continue
            fixed_weak += c*high_weak(x,p0,p1); others.append(F.lower())
            decay = [(-z/x).exp() for _,_,z in field]
            dk = c*F*sum(a*e for (_,a,_),e in zip(field,decay,strict=True))/pi
            dw = -c*F**3*sum(g*e for (g,_,_),e in zip(field,decay,strict=True))/(2*pi)
            row = [d*M[j,k].conjugate() for k,d in enumerate(delta)]
            for i,derivative in enumerate(dM):
                df = 2*sum((z*derivative[j,k]).real for k,z in enumerate(row))
                J[8,i] += dk*df; J[9,i] += dw*df
        return arb_mat([[v] for v in vals]), J, fixed_weak, min(others)

    point = [arb(0)]*ng
    for _ in range(5):
        F,J,_,_ = system(point,(t0+t1)/2)
        point = [(x-d).mid() for x,d in zip(point,J.solve(F).entries(),strict=True)]
    radius = arb('1e-7'); box=[x+arb(0,radius.upper()) for x in point]
    t=t0.union(t1); F,J,_,_=system(point,t); pre=J.inv()
    assert not pre.det().contains(0)
    _,wide,background,minimum = system(box,t,derivative_only=True)
    defect=arb_mat([[int(i==j) for j in range(ng)] for i in range(ng)])-pre*wide
    ratio=max((abs(v)/radius+sum(abs(x) for x in row)).upper()
              for v,row in zip((pre*F).entries(),defect.tolist(),strict=True))
    ends=[weights[15]*weak(z)+background for z in (t0,t1)]
    print('混合站、十面误差占半宽、其余振幅下界、两端弱差',sites[:len(seed[2])],
          ratio.str(9),minimum.str(9),*[v.str(9) for v in ends],flush=True)
    assert ratio<1 and ends[0]>0>ends[1]


def verify(root):
    shell=run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    tr=run_path(str(root/'ch39-adjacent-pvm-transport.py'))
    pv=run_path(str(root/'ch39-native-rh-ray-pvm.py'))
    gravity=run_path(str(root/'ch39-three-index-free-atoms.py'))
    high_weak=run_path(str(root/'ch39-full-neutral-free-cutoff.py'))['weak_column']
    nodes,pos=run_path(str(root/'ch39-left-cutoff-two-turns.py'))['verify'](root)
    assert all(pos[i] is pos[0] for i in range(19) if i not in ACTIVE)
    assert all(pos[i] is pos[17] for i in ACTIVE[:3])
    amps=[pos[i].sqrt() for i in (0,17,18)]; assert amps[0]<amps[1]<amps[2]
    n=CENTER-2;moments,heat,fields=gravity['broken_fields'](root,[n]);ctx.prec=320
    levels,source=tr['spectrum'](shell,6,n+1);_,target=tr['spectrum'](shell,6,n)
    assert all(x.overlaps(y) for x,y in zip(nodes,levels,strict=True))
    _,_,cross,_=next(v for v in pv['shell_weights'](root,n) if v[0]==6);ctx.prec=320
    weights=[z.real for z in cross[0]]; assert sum(weights).contains(1)
    floor=arb('91.1876')/arb('1.22089e19');four_site_bound(source,target,amps[1],floor)
    weak=pv['weak_series']();bounds=run_path(str(root/'ch39-outer-shared-distance-kernel.py'))['weak_bands']
    whole,large=bounds(weak,high_weak,floor**2,(arb('1.176397e14')/arb('1.22089e19'))**2)
    assert sum(c*(large if j==15 else whole) for j,c in enumerate(weights))<arb('-.01710')
    for seed in SEEDS:witness(seed,source,target,amps,weights,fields[0],moments,heat,weak,high_weak)


if __name__=='__main__':verify(Path(__file__).resolve().parent)
