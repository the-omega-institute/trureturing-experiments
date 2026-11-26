"""右外壳三荷幂次发生元的单射线障碍及双特征值谱宽必要界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


CENTER = 2_306_888_100_085
THRESHOLD = '6.54e-6'
SHIFT = '1.19308e-6'


def adj(x):
    return x.conjugate().transpose()


def prove_up(root, transport, shells, values):
    old = transport['spectrum'](shells,-4,CENTER+1)[1]
    new = transport['spectrum'](shells,-4,CENTER+2)[1]
    change = adj(old)*new
    h = adj(change)*acb_mat([[w*change[i,j] for j in range(19)]
                            for i,w in enumerate(values)])
    diag = arb_mat([[h[i,i].real] for i in range(19)])
    # p^* minors p 为二阶基本对称多项式。每个非对角项在二次型中出现两次。
    minors = arb_mat([[(diag[i,0]*diag[j,0]-abs(h[i,j])**2)/2 if i!=j else arb(0)
                       for j in range(19)] for i in range(19)])
    sites = [shells['soft']['pole'],*shells['boundary_sites']()]
    assert sites == [14,8,7,6]
    zero = [new[sites[0],j] for j in range(19)]
    one = [sum(new[k,j] for k in sites[1:])/arb(3).sqrt() for j in range(19)]
    pair = [u*v.conjugate() for u,v in zip(zero,one,strict=True)]
    coeff = arb_mat([[(abs(u)**2+abs(v)**2)/2 for u,v in zip(zero,one,strict=True)],
                     [v.real for v in pair], [v.imag for v in pair],
                     [(abs(u)**2-abs(v)**2)/2 for u,v in zip(zero,one,strict=True)]])
    target = arb(THRESHOLD)
    q = coeff*minors*coeff.transpose()
    mean = coeff*diag
    for i in range(4):
        q[0,i] -= target*mean[i,0]/2
        q[i,0] -= target*mean[i,0]/2
        q[i,i] += (-1 if i==0 else 1)*arb(SHIFT)
    leading = [arb_mat([[q[i,j] for j in range(k)] for i in range(k)]).det() for k in range(1,5)]
    assert all(v>0 for v in leading)
    # 标量输运时同一平凡二维块的最小值是 mean_0 - |mean_xyz|。
    raw = [(mean[0,0]-sum(mean[i,0]**2 for i in range(1,4)).sqrt()).lower()]
    omega = acb(arb('-0.5'),arb(3).sqrt()/2)
    for z in (omega,omega.conjugate()):
        v = [sum(z**k*new[i,j] for k,i in enumerate(sites[1:]))/arb(3).sqrt() for j in range(19)]
        p = arb_mat([[abs(a)**2] for a in v])
        first = (p.transpose()*diag)[0,0]
        second = (p.transpose()*minors*p)[0,0]
        assert second > target*first
        raw.append(first.lower())
    assert min(raw) > arb('.0019196') and max(v.upper() for v in values) < arb('.94756')
    # 对任意秩，用两个本征相位的复平面弦中点作标量比较。
    s = arb('.042')
    assert (1-s**2).sqrt()*arb('.0019196').sqrt()-s*arb('.94756').sqrt() > arb('6.29e-6').sqrt()
    print('单射线上型下界、四个顺序主子式、无转动时四站下界、任意秩必要的 |sin θ|',
          THRESHOLD,[m.str(10) for m in leading],min(raw).str(12),s.str(5),flush=True)


def prove_weak(transport, shells, training):
    old_levels, old = transport['spectrum'](shells,6,CENTER+1)
    new = transport['spectrum'](shells,6,CENTER+2)[1]
    weights = transport['training_values'](old_levels,training)
    selector = adj(old)*acb_mat([[new[i,15]] for i in range(19)])
    cutoff = 16*arb('91.1876')/arb('1.22089e19')
    high = [j for j,w in enumerate(weights) if w>cutoff]
    assert high == [17,18] and all(w<cutoff for w in weights[:17])
    m = min(w.lower() for w in weights)
    h = min(weights[j].lower() for j in high)
    leakage = ((cutoff-m)/(h-m)).sqrt()
    initial = sum(abs(selector[j,0])**2 for j in high).sqrt()
    sine = initial*(1-leakage**2).sqrt()-(1-initial**2).sqrt()*leakage
    assert sine > arb('.23454770')
    print('第十六层源高投影、D=0允许的高投影、|sin(3θ/2)| 严格下界',
          initial.str(13),leakage.str(9),sine.lower().str(13),flush=True)


def verify(root):
    training,_ = run_path(str(root/'ch39-three-shell-two-space-readout.py'))['verify'](root,return_data=True)
    shells = run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    transport = run_path(str(root/'ch39-adjacent-pvm-transport.py'))
    ctx.prec = 384
    levels = transport['spectrum'](shells,-4,CENTER+1)[0]
    up = transport['training_values'](levels,training)
    prove_up(root,transport,shells,up)
    prove_weak(transport,shells,training)


if __name__=='__main__':
    verify(Path(__file__).resolve().parent)
