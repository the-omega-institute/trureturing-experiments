"""同一枚实际零点的重偶倍率、钟弦分母与两圈规范及配平曲率的共同窗。"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


def verify(root):
    physical = run_path(str(root / 'ch39-first-three-mode-observables.py'))
    model, initial = physical['model_inputs']()
    ctx.prec, ctx.cap = 256, 3
    joint = physical['joint']
    vacuum = joint['vacuum']
    chiral = joint['chiral']
    height = acb.zeta_zero(2_306_888_100_085).imag
    rows = chiral['source']['spectral_supplier'](model)(height)
    layers = {q:chiral['ordered']([v for v, charge in rows if charge == q]) for q in (0, *model['charges'])}
    assert all(layers[q][2]<9<10<layers[q][3] for q in model['charges'])
    k0, d0 = arb(32), arb(724)
    mp, reduced, vev = map(arb, ('1.22089e19','2.435e18','246.22'))
    coefficient = (9*arb(3).log()-12*arb(2).log())/2
    target = arb('1e-122')
    a0 = (vev**2/(mp*(k0*reduced/d0**2)*layers[2][15]))**2
    main = coefficient*a0**2/(arb.pi()*target)
    def scaled_rows(k,d):
        # fields 以 32 和 724 为单位；输入谱重标恰好替换其中每一行实际质量。
        bands = {q:[v*(d0/d)**2*(k/k0 if q%2==0 else 1) for v in band]
                 for q,band in layers.items()}
        original = vacuum['minimal']['fields'](bands)
        # 与原处方相同，每个有质量粒子分别配平。相遇的行可直接相加。
        residual = [(weight+(arb(-g)/6 if g<0 else arb(g)/12), z) for g,weight,z in original]
        a = a0*(k0/k)**2*(d/d0)**4
        return residual+vacuum['low_sector'](4,a), a

    kbox = arb('31.98').union(arb('32.02'))
    dbox = arb('723.8').union(arb('724.2'))
    fields,a = scaled_rows(kbox,dbox)
    positive = sum(w for w,z in fields if w>0)
    negative = sum(w for w,z in fields if w<0)
    # 对任意正 k,d，权与参数无关。牛顿标定处 E2 <= 1 给全域粗下界。
    xmin = arb('.006')
    assert 2*arb.pi()/positive > xmin
    assert coefficient*(2*a0)**2/(arb.pi()*target)*(1-2*a0/(xmin*coefficient)) > arb('1.002')
    assert main/4 < arb('.998')
    # 曲率的原积分对 a,x 均递增，故过窗必有 a0/2<a<2a0。
    loss = 2*a0/(xmin*coefficient)
    assert loss < arb('1e-57')
    low, high = arb('.998'),arb('1.002')
    # I(0)-u <= I(u) <= I(0)，每个端点夹在下式的相邻两端内。
    effective = main*(1-arb(0).union(loss))
    kends = [k0*(effective/bound).root(4) for bound in (high,low)]
    dends = [d0*(bound/effective).root(8) for bound in (low,high)]
    assert all(kbox.contains(z) for z in kends) and all(dbox.contains(z) for z in dends)
    assert kbox.lower()>31 and kbox.upper()<33 and dbox.lower()>723 and dbox.upper()<725
    # 有界参数域上的牛顿根存在，且在全域下界 xmin 以上惟一。
    assert vacuum['newton'](fields,arb('.03')) < 1 < vacuum['newton'](fields,arb('.033'))
    derivative_floor=sum(w*(-z.upper()/xmin).exp() if w>0 else w for w,z in fields)
    assert derivative_floor>100, derivative_floor
    # 由 (x E2(z/x))'=exp(-z/x)、u E1(u)<=1 隐式求导。
    cutoff_slope = 4*(positive-negative)/derivative_floor
    assert cutoff_slope<25
    # |I'(u)|<=1，质量平方对 log k / log d 的斜率依次为 -2 / 4（低层）。
    slope_error=(4+cutoff_slope)*a/(arb('.03')*coefficient-a)
    assert slope_error < arb('1e-56')
    # 原跑动中只包含奇荷开启的阈值，偶荷脱耦对变化后的 k 另作核验。
    even = [*layers[0],*(v for q in (2,4,6) for v in layers[q][3:])]
    mu,alpha=arb(0),arb(0)
    slices=8
    for i in range(slices):
        left=arb('723.8')+arb('.4')*i/slices
        right=arb('723.8')+arb('.4')*(i+1)/slices
        d=left.union(right)
        with redirect_stdout(StringIO()):
            energy,_,strong=joint['gauge']['crossing_enclosure'](model,initial,layers,False,d,'29.04','.1184',recenter=True)
        assert arb('.117')<strong<arb('.119') and all(kbox*v*reduced/d**2>energy for v in even)
        mu,alpha=(energy,strong) if i==0 else (mu.union(energy),alpha.union(strong))
    # 两个端面外的曲率已越界，连续、严格单调给出各自的整个单参数区间。
    assert effective*(k0/kbox.lower())**4 > high and effective*(k0/kbox.upper())**4 < low
    assert effective*(dbox.lower()/d0)**8 < low and effective*(dbox.upper()/d0)**8 > high
    # 二者同调可同时越过两条单参数窗口的上端，而保持跷跷板质量。
    mixed_d = arb('724.19')
    mixed_k = k0*(mixed_d/d0)**2
    assert mixed_d > dends[1] and mixed_k > kends[1]
    assert dbox.contains(mixed_d) and kbox.contains(mixed_k) and low<effective<high
    strong_edges=[]
    for d in (arb('723.82'),arb('724.18')):
        with redirect_stdout(StringIO()):
            strong_edges.append(joint['gauge']['crossing_enclosure'](model,initial,layers,False,d,'29.04','.1184',recenter=True)[2])
    secant=(strong_edges[1]-strong_edges[0])/arb('.36')
    assert secant<0
    print('保持重偶质量的双参数例与分母两点的强耦合割线斜率',mixed_d.str(8),mixed_k.str(16),secant.lower().str(9),secant.upper().str(9))
    print('实际零点谱的第十六个荷二质量、名义曲率主项', layers[2][15].str(18), main.str(18))
    print('各端点：固定分母时的重偶倍率，固定倍率时的分母', *(v.str(19) for v in [*kends,*dends]))
    print('整域外包的强耦合、交点能标 / GeV', alpha.str(10),mu.str(10))
    print('曲率对倍率、分母的对数导数分别距 -4、8 小于', slope_error.abs_upper().str(8))
    print('阈值处最轻偶质量与最大交点之比下界', (kbox.lower()*min(v.lower() for v in even)*reduced/dbox.upper()**2/mu.upper()).str(9))


if __name__=='__main__':
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
