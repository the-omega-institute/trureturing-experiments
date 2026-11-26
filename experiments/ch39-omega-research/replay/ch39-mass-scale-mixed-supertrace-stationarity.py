"""两个质量标尺 (k,d) 上混合超迹泛函的驻点：可分离约化、四维截止族驻点个数与分支的完全分类与两幂族。

说明见“常数锚混合超迹泛函的驻点分母由锚决定而四维截止族须调准一个无量纲系数才能进入两质量标尺窗.md”。

沿用 ch39-mass-scale-spectral-stationarity.py 的场表与记号：98 行有质量重场，奇荷 (1,3) 质量
Mbar*h/d^2，偶荷 (0,2,4,6) 质量 Mbar*k*h/d^2，低能四层不计入。锚标度 Lambda 为与 (k,d) 无关的常数。
记 y_o=(Mbar/Lambda)^2/d^4，y_e=k^2 y_o，则 x_r=m_r^2/Lambda^2=y_s h_r^2，k=sqrt(y_e/y_o)。

一、四维截止族  F/Lambda^4 = sum_r g_r [x_r^2 (log x_r - c) + kappa x_r]，不变量 sigma=kappa e^{1/2-c}。
    每个扇区的驻定方程化为 w e^w = -sigma A_s/2，A_s=(S2/S4) exp(L4/S4)，L4=sum g h^4 log h^2。
二、两幂族  F = STr M^{2p} + lam Lambda^{2p-2q} STr M^{2q}（整数 1<=q<p）：k_* 与 lam 无关。
三、奇偶分离的减法常数与锚：k_* = k_0 exp((c_e-c_o)/2) = k_0 Lambda_e/Lambda_o。
"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx


def inside(x, lo, hi):
    assert arb(lo) < x < arb(hi), (x, lo, hi)
    return x


def verify(root):
    prev = run_path(str(root/'ch39-mass-scale-spectral-stationarity.py'))
    physical = run_path(str(root/'ch39-first-three-mode-observables.py'))
    model, _ = physical['model_inputs']()
    ctx.prec, ctx.cap = 192, 3
    joint = physical['joint']
    chiral, vacuum = joint['chiral'], joint['vacuum']
    layers = prev['layers_at'](chiral, model, prev['FLAGSHIP'])
    rows = prev['field_rows'](vacuum, layers)          # 逐行核对原 minimal['fields']
    odd = [(arb(g), h) for g, _, e, h in rows if not e]
    even = [(arb(g), h) for g, _, e, h in rows if e]
    assert len(odd) == 32 and len(even) == 66 and all(g < 0 for g, _ in odd)

    def S(sec, p):
        return sum(g*h**(2*p) for g, h in sec)

    def L(sec, p):
        return sum(g*h**(2*p)*(h*h).log() for g, h in sec)

    sec = {'o': odd, 'e': even}
    S2 = {s: S(v, 1) for s, v in sec.items()}
    S4 = {s: S(v, 2) for s, v in sec.items()}
    off = {s: L(v, 2)/S4[s] for s, v in sec.items()}   # = 2<log h>_{s, g h^4}
    A = {s: S2[s]/S4[s]*off[s].exp() for s in sec}
    assert S2['o'] < 0 and S4['o'] < 0 and S2['e'] > 0 and S4['e'] > 0
    inside(A['o'], '1.45097670', '1.45097671')
    inside(A['e'], '1.1008010', '1.1008012')
    assert A['o'] > 1 and A['o'] > A['e'] > 0             # A_o>=1 是 Jensen 不等式（奇荷权同号）
    sigc = {s: 2/(arb(1).exp()*A[s]) for s in sec}       # 扇区有零点当且仅当 sigma<=sigc_s
    inside(sigc['o'], '0.50707835', '0.50707836')
    inside(sigc['e'], '0.6683849', '0.6683851')
    k0 = ((off['o']-off['e'])/2).exp()
    inside(k0, '2.4134605455', '2.4134605465')            # 与上一证书的四维 CW 倍率一致
    print('A_o, A_e', A['o'].str(12), A['e'].str(12), '临界 sigma', sigc['o'].str(12), sigc['e'].str(12))

    # 两窗与 k>=1 给出的必要倍率区间（上一证书的外端），以及 d=724 处的曲率必要带外端。
    klo, khi = arb('29.9841147'), arb('35.0010644')
    blo, bhi = arb('31.98396690'), arb('32.01596692')

    def logk(sigma, bo, be):
        wo = (-sigma*A['o']/2).lambertw(branch=bo)
        we = (-sigma*A['e']/2).lambertw(branch=be)
        assert wo.is_finite() and we.is_finite()
        return (we-wo-off['e']+off['o'])/2, wo, we

    # 一、四维截止族的分支包络（书面证明给出不等式，这里给出数值端点）。
    kpoly = (S2['e']*S4['o']/(S4['e']*S2['o'])).sqrt()
    assert (kpoly-k0*(A['e']/A['o']).sqrt()).contains(0)
    inside(kpoly, '2.1021534145', '2.1021534146')        # sigma<0：k 在 (kpoly, k0)
    top00 = k0*arb(0.5).exp()
    inside(top00, '3.979123', '3.979125')          # (W0,W0)：k 在 (k0, k0 e^{1/2})
    assert top00 < klo
    # sigma=sigc_o 处（奇向重根 w_o=-1）的端点值。
    ze = -(A['e']/A['o'])/arb(1).exp()
    kc_m1_0 = inside(((ze.lambertw()+1-off['e']+off['o'])/2).exp(), '3.212017', '3.212019')
    kc_0_m1 = inside(((ze.lambertw(branch=-1)+1-off['e']+off['o'])/2).exp(), '1.510152', '1.510154')
    # 分支 (W0,W-1) 与 (W-1,W-1) 的 k<=k0；分支 (W-1,W0) 的 k 从 +inf 单调降到 kc_m1_0。
    # 抽样核对单调性与分支不等式（全称由书面证明给出）。
    grid = [sigc['o']*arb(i)/64 for i in range(1, 64)]
    prev_k = None
    for s in grid:
        a00, a01, a10, a11 = (logk(s, bo, be)[0] for bo, be in ((0, 0), (0, -1), (-1, 0), (-1, -1)))
        assert k0.log() < a00 < top00.log() and a01 < k0.log() and a11 < k0.log()
        assert kc_m1_0.log() < a10
        if prev_k is not None:
            assert a10 < prev_k
        prev_k = a10
    for s in (arb(-1)/1000, arb(-1), arb(-1000), -arb(10)**6):
        a = logk(s, 0, 0)[0]
        assert kpoly.log() < a < k0.log()
    print('sigma<0 的 k 包络', kpoly.str(11), k0.str(11), '(W0,W0) 上界', top00.str(10),
          'sigma=sigc_o 处分支端点', kc_m1_0.str(11), kc_0_m1.str(11))

    # 分支 (W-1 奇, W0 偶)：log k 关于 sigma 严格递减，用二分求目标倍率的原像括号。
    def preimage(T):
        a, b = arb('0.02'), arb('0.06')
        assert logk(a, -1, 0)[0] > T.log() > logk(b, -1, 0)[0]
        for _ in range(40):                                   # 用中点值定位，最后再严格认证括号
            m = (a+b)/2
            if logk(m, -1, 0)[0].mid() > T.log().mid():
                a = m
            else:
                b = m
        step = arb('1e-11')
        while not (logk(a, -1, 0)[0] > T.log() > logk(b, -1, 0)[0]):
            a, b, step = a-step, b+step, 2*step
            assert step < arb('1e-6')
        ka, kb = logk(a, -1, 0)[0], logk(b, -1, 0)[0]
        assert ka > T.log() > kb                              # 严格括号（单调性见书面证明）
        return a, b
    expect = {'khi': ('0.0345056616', '0.0345056620'), 'bhi': ('0.0397731745', '0.0397731749'),
              '32': ('0.0398046438', '0.0398046442'), 'blo': ('0.0398362825', '0.0398362829'),
              'klo': ('0.0441123193', '0.0441123197')}
    pre = {}
    for name, T in (('khi', khi), ('bhi', bhi), ('32', arb(32)), ('blo', blo), ('klo', klo)):
        a, b = preimage(T)
        lo, hi = expect[name]
        assert arb(lo) < a and b < arb(hi)
        pre[name] = (a, b)
    # 窗：sigma 在 (sigma(khi), sigma(klo)) 内是 k 落入必要区间的充要条件；曲率带对应更窄的区间。
    win_rel = inside(pre['klo'][1]/pre['khi'][0], '1.2783', '1.2785')
    band_rel = inside(pre['blo'][1]/pre['bhi'][0]-1, '0.0015866', '0.0015868')
    s32 = pre['32'][0].union(pre['32'][1])
    print('k 落入 (29.9841147,35.0010644) 的 sigma 区间', pre['khi'][0].str(10), pre['klo'][1].str(10))
    print('k 落入曲率带的 sigma 区间', pre['bhi'][0].str(11), pre['blo'][1].str(11), '相对宽度', band_rel.str(5))
    print('k_*=32 的 sigma', s32.str(11))
    # 仅数值（依赖上一证书的点证书）：滑行点 d=716 与 d=743 的倍率在窗内，对应的 sigma 原像。
    slide = {}
    for name, T in (('743', arb('33.7015964104')), ('716', arb('31.2967247642'))):
        slide[name] = preimage(T)
    inside(slide['743'][0], '0.03665655', '0.03665665')
    inside(slide['716'][1], '0.04122995', '0.04123025')
    slide_rel = inside(slide['716'][1]/slide['743'][0]-1, '0.1247', '0.1249')
    print('滑行点 716、743 对应的 sigma', slide['743'][0].str(9), slide['716'][1].str(9), '相对宽', slide_rel.str(4))
    # sigma=1/25 在必要区间内但不在曲率带内（仅作对照）。
    k25 = logk(arb(1)/25, -1, 0)[0].exp()
    inside(k25, '31.90131', '31.90133')
    assert klo < k25 < blo

    # 在 sigma_32、c=3/2 处独立核对驻点：直接按行求一阶导数含零，二阶导数两向为正（局部极小）。
    c = arb(3)/2
    _, wo, we = logk(s32, -1, 0)
    ly = {'o': wo+c-arb(1)/2-off['o'], 'e': we+c-arb(1)/2-off['e']}
    kappa32 = inside(s32*(c-arb(1)/2).exp(), '0.10820024', '0.10820025')
    for s in ('o', 'e'):
        y = ly[s].exp()
        xs = [(g, y*h*h) for g, h in sec[s]]
        first = sum(g*(x*x*(2*(x.log()-c)+1)+kappa32*x) for g, x in xs)
        second = sum(g*(x*x*(4*(x.log()-c)+4)+kappa32*x) for g, x in xs)
        size = sum(abs(g)*(x*x*abs(x.log()-c)+x) for g, x in xs)
        assert first.contains(0) and abs(first) < size*arb('1e-6')
        assert second > 0
    kcheck = ((ly['e']-ly['o'])/2).exp()
    assert kcheck.contains(32)
    mbar = arb('2.435e18')
    lam724 = inside(mbar/(arb(724)**2*(ly['o']/2).exp()), '5.6292538e15', '5.6292542e15')
    xmax = {}
    for t in ('o', 'e'):
        vals = [ly[t].exp()*h*h for _, h in sec[t]]
        m = vals[0]
        for v in vals[1:]:
            m = m.max(v)                                     # 逐个归约，不比较可能重叠的球
        xmax[t] = m
    xo_max = inside(xmax['o'], '0.0198874', '0.0198876')
    xe_max = inside(xmax['e'], '5.1316', '5.1317')
    lam_ratio = inside(bhi/blo, '1.0010004', '1.0010006')
    print('c=3/2 时 kappa_32', kappa32.str(10), '使 d_*=724 的 Lambda/GeV', lam724.str(8),
          '奇/偶最大 x', xo_max.str(6), xe_max.str(6), 'Lambda 可容许比', lam_ratio.str(8))
    # 规范截止的截断（c=1/2, kappa=2 即 sigma=2）落在无驻点区。
    assert arb(2) > sigc['e'] > sigc['o']
    # 固有时截止的大 Lambda 截断：c=3/2-gamma、kappa=2，sigma=2 e^{gamma-1}，同样落在无驻点区。
    sig_pt = inside(2*(arb.const_euler()-1).exp(), '1.31043', '1.31045')
    assert sig_pt > sigc['e']
    # 按场型取不同减法常数（MS-bar 的 3/2 与 5/6，跨度 2/3）：沿用上一证书的跨度界
    # |log k* - log k0| <= Delta (R+1)/4，R = sum_e |g| h^4 / S4_e。
    R = sum(abs(g)*h**4 for g, h in even)/S4['e']
    inside(R, '10.212445', '10.212447')
    span = arb(3)/2-arb(5)/6
    assert span < arb('0.8988604')
    k_span = inside(k0*(span*(R+1)/4).exp(), '15.6395', '15.6396')
    assert k_span < klo
    print('固有时截断 sigma', sig_pt.str(8), '跨度 2/3 时 k* 上界', k_span.str(6))

    # 二、两幂族：整数 1<=q<p<=30 逐对计算；p>30 由上一证书的尾部包络与凸组合处理。
    Sp = {p: (S(odd, p), S(even, p)) for p in range(1, 31)}
    assert all(Sp[p][0] < 0 for p in Sp)
    assert all((Sp[p][1] > 0) if p <= 5 else (Sp[p][1] < 0) for p in Sp)
    table = {}
    for p in range(2, 31):
        for q in range(1, p):
            ro, re = Sp[q][0]/Sp[p][0], Sp[q][1]/Sp[p][1]
            assert ro > 0
            if (p <= 5) != (q <= 5):
                assert re < 0                          # 跨越偶荷分母变号：无驻点
                continue
            assert re > 0
            table[(p, q)] = (re/ro).log()/(2*(p-q))
    assert len(table) == 10+300
    best = (5, 4)
    assert all(table[t] < table[best] for t in table if t != best)
    kbest = inside(table[best].exp(), '3.3260466', '3.3260467')
    kmin = None
    for t in table.values():
        kmin = t if kmin is None else kmin.min(t)
    inside(kmin.exp(), '0.9580836', '0.9580838')
    assert (kmin-table[(7, 6)]).contains(0) and all(table[t] > table[(7, 6)] for t in table if t != (7, 6))
    print('两幂族整数对的最大倍率 (p,q)=(5,4)', kbest.str(11), '最小', kmin.exp().str(8),
          '可取驻点的对数', len(table))
    tail_hi = arb('2.027536')                            # 上一证书：p>=30 实数尾部 k_CW<2.027536
    assert kbest > tail_hi and kbest < klo
    # 由 log k(p,q) = [(30-q)log k(30,q) + (p-30) log k(p,30)]/(p-q)，p>30>=q>=6 的值不超过 max(k(30,q), tail_hi)。
    assert all(table[(30, q)] < kbest.log() for q in range(6, 30))
    assert (table[(2, 1)].exp()-kpoly).contains(0)           # (2,1) 对即 sigma->-inf 的极限

    # 三、奇偶分离的减法常数：Delta c = c_e-c_o 与 log k 一一对应。
    dc_lo = inside(2*(klo/k0).log(), '5.03921', '5.03923')
    dc_hi = inside(2*(khi/k0).log(), '5.348633', '5.348635')
    dc_32 = inside(2*(arb(32)/k0).log(), '5.169348', '5.169350')
    ratio_32 = inside(arb(32)/k0, '13.25896', '13.25898')
    dc_band = inside(2*(bhi/blo).log(), '0.0019999', '0.0020001')
    print('分离常数：窗要求 c_e-c_o 在', dc_lo.str(7), dc_hi.str(7), 'k=32 时', dc_32.str(7),
          '曲率带宽', dc_band.str(5), '分离锚之比', ratio_32.str(8))


if __name__ == '__main__':
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
