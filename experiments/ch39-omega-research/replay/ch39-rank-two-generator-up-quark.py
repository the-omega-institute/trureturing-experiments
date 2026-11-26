"""右侧相邻壳用同一秩二发生元的上型幂次精确实现全部三个上型 Yukawa。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


N = 2_306_888_100_085
CENTERS = ('.00000615', '.00311', '.8616')
WINDOWS = (('.00000601','.00000629'),('.00306','.00316'),('.8573','.8659'))
REAL = ((-2655,-1965),(823,-4327),(-420,-488),(3764,-472),(2828,3841),(-4007,3953),(237,54))
IMAG = ((1700,4346),(1590,3257),(-1374,-1143),(-3100,-2147),(-3688,-43),(4826,-3329),(-15,281))


def adj(a):
    return a.conjugate().transpose()


def eye(n):
    return acb_mat([[int(i==j) for j in range(n)] for i in range(n)])


def projector(orbits):
    basis = acb_mat([[1/arb(len(g)).sqrt() if i in g else 0 for g in orbits] for i in range(19)])
    seed = acb_mat([[acb(x,y) for x,y in zip(a,b,strict=True)] for a,b in zip(REAL,IMAG,strict=True)])
    assert (adj(seed)*seed).det().real > 0
    z = basis*seed
    p = z*(adj(z)*z).inv()*adj(z)
    assert (p*p-p).contains(acb_mat(19,19)) and p.trace().contains(2)
    return p


def read(source, target, values, u):
    overlap = adj(target)*u*source
    return [sum(abs(overlap[j,i])**2*w for i,w in enumerate(values)) for j in range(19)]


def compression(raw, masses):
    values,eigen = raw.eig(right=True,algorithm='rump')
    order = sorted(range(4),key=lambda i:values[i].real.mid())
    lam = [values[i].real for i in order]
    assert all(lam[j] < arb(lo) < arb(hi) < lam[j+1] for j,(lo,hi) in enumerate(WINDOWS))
    columns = [acb_mat([[eigen[i,j]] for i in range(4)]) for j in order]
    columns = [v*(1/(adj(v)*v)[0,0].real.sqrt()) for v in columns]
    q = acb_mat([[v[i,0] for v in columns] for i in range(4)])
    # 在四阶本征基中，截取压缩特征多项式/原特征多项式的正留数。
    weights = [prod(a-m for m in masses)/prod(a-b for k,b in enumerate(lam) if k!=i)
               for i,a in enumerate(lam)]
    assert min(weights)>0 and sum(weights).contains(1)
    normal = [x.sqrt() for x in weights]
    vectors = []
    for m in masses:
        v = acb_mat([[n/(a-m)] for n,a in zip(normal,lam,strict=True)])
        vectors.append(q*v*(1/(adj(v)*v)[0,0].real.sqrt()))
    frame = acb_mat([[v[i,0] for v in vectors] for i in range(4)])
    assert (adj(frame)*frame).contains(eye(3))
    diag = acb_mat([[masses[i] if i==j else 0 for j in range(3)] for i in range(3)])
    diff = adj(frame)*raw*frame-diag
    residual = max(abs(a).upper() for a in diff.entries())
    assert diff.contains(acb_mat(3,3)) and residual < arb('1e-16'), (residual, diff.contains(acb_mat(3,3)))
    print('四站上型振幅本征值',[v.str(16) for v in lam],flush=True)
    print('四站法向在有序本征基中的概率',[v.str(12) for v in weights],flush=True)
    return frame


def verify(root):
    training,_ = run_path(str(root/'ch39-three-shell-two-space-readout.py'))['verify'](root,return_data=True)
    shells = run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    transport = run_path(str(root/'ch39-adjacent-pvm-transport.py'))
    orbits = run_path(str(root/'ch39-up-quark-arbitrary-site-phase-floor.py'))['ORBITS']
    ctx.prec = 384
    p = projector(orbits)
    x = eye(19)+(acb(0,1)-1)*p
    up = eye(19)-2*p
    gravity = eye(19)+(acb(0,-1)-1)*p
    assert (adj(x)*x).contains(eye(19)) and (x*x*up).contains(eye(19))
    assert (x*x*x-gravity).contains(acb_mat(19,19))
    assert (1/arb(2).sqrt()) > arb('.23454770')
    for charge,u in ((-4,up),(6,gravity)):
        levels,source = transport['spectrum'](shells,charge,N+1)
        _,target = transport['spectrum'](shells,charge,N+2)
        amplitudes = read(source,target,transport['training_values'](levels,training),u)
        if charge==-4:
            sites = (shells['soft']['pole'],*shells['boundary_sites']())
            c = acb_mat([[target[i,j] for j in range(19)] for i in sites])
            raw = acb_mat([[c[i,j]*amplitudes[j] for j in range(19)] for i in range(4)])*adj(c)
            frame = compression(raw,list(map(arb,CENTERS)))
        else:
            cutoff = 16*arb('91.1876')/arb('1.22089e19')
            assert amplitudes[15]>cutoff
            print('同一生成元荷六第十六层振幅及零弱差要求的上限',
                  amplitudes[15].str(13),cutoff.str(10),flush=True)
    return x,frame


if __name__=='__main__':
    verify(Path(__file__).resolve().parent)
