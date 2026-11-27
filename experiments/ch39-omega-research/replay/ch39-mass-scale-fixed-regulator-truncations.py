"""两个质量标尺 (k,d) 上由正则方案固定系数的超迹泛函：正 Schwinger 权类的截断下界与四维加六维 CW 混合族的完全分类。

说明见“取值不超过1的正Schwinger权正则的对数阶截断给出不变量不低于一点三一零四而四维加六维CW混合在两种常数取法下须把系数调进两个窄区间才能进入两质量标尺窗.md”。

沿用 ch39-mass-scale-mixed-supertrace-stationarity.py 的场表与记号：98 行有质量重场，奇荷 (1,3) 质量
Mbar*h/d^2，偶荷 (0,2,4,6) 质量 Mbar*k*h/d^2，低能四层不计入；锚标度 Lambda 为与 (k,d) 无关的常数，
y_o=(Mbar/Lambda)^2/d^4，y_e=k^2 y_o，x_r=y_s h_r^2，k=sqrt(y_e/y_o)。

一、正 Schwinger 权类  f'(x)=int_0^inf s^{-2} R(s) e^{-sx} ds，R>=0。对数阶截断 f'~M0+x log x-C x 即四维截止族
    成员，sigma=2 M0 e^{-C}。浴缸不等式给出 sigma>=2 e^{gamma-1-E}，E=int (R-1)_+ ds/s；R<=1 时 sigma>=2e^{gamma-1}。
二、四维加六维 CW 混合  x^2(log x-c4)+lam x^3(log x-c6)：不变量 tau=lam e^{c4-1/2}、delta=c6-c4；每个扇区的驻定
    方程写成 tau=T_s(t)=-K_s e^{-t}(t-tz_s)/(t-tp_s)，t=log y（以 e^{c4-1/2} 为单位），单调分段可显式写出。
    对 delta=0 与 delta=1/3，用单调分段的原像集合逐格认证 tau 轴上进入两窗必要区间的全部可能位置。
"""
from pathlib import Path
from runpy import run_path
import random
import sys

from flint import arb, ctx


def inside(x, lo, hi):
    assert arb(lo) < x < arb(hi), (x, lo, hi)
    return x


def sector_sums(root):
    prev = run_path(str(root/'ch39-mass-scale-spectral-stationarity.py'))
    physical = run_path(str(root/'ch39-first-three-mode-observables.py'))
    model, _ = physical['model_inputs']()
    ctx.prec, ctx.cap = 256, 3
    joint = physical['joint']
    layers = prev['layers_at'](joint['chiral'], model, prev['FLAGSHIP'])
    rows = prev['field_rows'](joint['vacuum'], layers)       # 逐行核对原 minimal['fields']
    sec = {'o': [(arb(g), h) for g, _, e, h in rows if not e], 'e': [(arb(g), h) for g, _, e, h in rows if e]}
    assert len(sec['o']) == 32 and len(sec['e']) == 66 and all(g < 0 for g, _ in sec['o'])
    S = {(s, p): sum(g*h**(2*p) for g, h in v) for s, v in sec.items() for p in (1, 2, 3)}
    L = {(s, p): sum(g*h**(2*p)*(h*h).log() for g, h in v) for s, v in sec.items() for p in (2, 3)}
    hmax2 = {}
    for k_, v in sec.items():
        m = None
        for _, h in v:
            m = h*h if m is None else m.max(h*h)              # 逐个归约，不比较可能重叠的球
        hmax2[k_] = m
    return S, L, hmax2


# ---------------------------------------------------------------- 一、正 Schwinger 权类
def step_sigma(knots, vals):
    """R 在 [knots[i],knots[i+1]) 上取 vals[i]，knots[0] 之前为 0，knots[-1] 之后为 1。返回 (sigma, E, M0)。"""
    M0 = sum(v*(1/a-1/b) for v, a, b in zip(vals, knots[:-1], knots[1:])) + 1/knots[-1]
    C = 1-arb.const_euler() + sum(v*(b/a).log() for v, a, b in zip(vals, knots[:-1], knots[1:])) - knots[-1].log()
    sstar = 1/M0
    E = arb(0)
    for v, a, b in zip(vals, knots[:-1], knots[1:]):
        if v > 1:
            lo = a.max(sstar) if sstar < b else None     # (R-1)_+ (1-s*/s)/s 在 [max(a,s*),b] 上的积分
            if lo is not None:
                E += (v-1)*((b/lo).log() + sstar*(1/b-1/lo))
    return 2*M0*(-C).exp(), E, M0


def schwinger_part(sig_o, sig_e, sig_win_hi):
    g = arb.const_euler()
    floor = inside(2*(g-1).exp(), '1.3104398516', '1.3104398517')
    assert floor > sig_e > sig_o                            # 截断成员一律落在两扇区都无零点的区
    # 闭式成员：R 为正随机变量 S 的分布函数时 sigma = 2 e^{gamma-1} E[1/S] exp(E log S)。
    hard = 2*(g-1).exp()*(1-g).exp()                        # 硬动量截止：S~Gamma(2)，E[1/S]=1，E log S=1-gamma
    assert (hard-2).contains(0)
    heat = inside(2*g.exp(), '3.5621', '3.5622')            # 热核动量截止 rho=e^{-t}：S 密度 2s/(1+s)^3
    pv = inside(8*arb(2).log()/arb(1).exp(), '2.039956', '2.039957')  # 最少 PV（质量平方 Lambda^2, 2Lambda^2）
    assert hard > floor and heat > floor and pv > floor      # 逐个比较
    # 独立数值核对上述闭式（中点积分，仅作核对）：M0=int s^{-2}R，C=1-gamma+int (R-1_{s>1})/s。
    import mpmath as mp
    mp.mp.dps = 30
    for name, R, want in (('hard', lambda s: mp.gammainc(2, 0, s, regularized=True), 2), ('heat', lambda s: s*s/(1+s)**2, 2*mp.e**mp.euler),
                          ('pv', lambda s: mp.expm1(-s)**2, 8*mp.log(2)/mp.e)):
        M0 = mp.quad(lambda s: R(s)/s**2, [0, 1, mp.inf])
        C = 1-mp.euler + mp.quad(lambda s: R(s)/s, [0, 1]) + mp.quad(lambda s: (R(s)-1)/s, [1, mp.inf])
        assert abs(2*M0*mp.e**(-C)-want) < mp.mpf(10)**-20, name
        # 推论（说明中定理 2 第 2 款）：R<=1 时截断 T=M0+x log x-Cx 处处不低于精确 f'；抽样核对
        for x in (mp.mpf('0.01'), mp.mpf('0.3'), mp.mpf(1), mp.mpf(5)):
            fp = mp.quad(lambda s: R(s)/s**2*mp.e**(-s*x), [0, 1, mp.inf])
            assert 0 < fp < M0 + x*mp.log(x) - C*x, (name, x)
    # 浴缸不等式的独立核对：随机阶梯 R（精确 Arb 算术）满足 sigma >= 2e^{gamma-1-E}；R<=1 时 sigma >= floor。
    rnd = random.Random(20261127)
    used = 0
    for trial in range(400):
        n = rnd.randint(1, 6)
        knots = sorted(arb(rnd.randint(1, 4000))/100 for _ in range(n+1))
        if any(not (b > a) for a, b in zip(knots[:-1], knots[1:])):
            continue
        cap = 1 if trial % 2 == 0 else 4
        vals = [arb(rnd.randint(0, 100*cap))/100 for _ in range(n)]
        sig, E, _ = step_sigma(knots, vals)
        used += 1
        # 定理给出的是 >=；等号在 R=1_[s*,inf) 时取到（随机样本可以恰好是它），故允许包围含零
        assert sig > 2*(g-1-E).exp() or (sig-2*(g-1-E).exp()).contains(0)
        if cap == 1:
            assert sig > floor or (sig-floor).contains(0)
    assert used == 400, used
    print('随机阶梯样本数', used)
    # 进窗的必要超调量：sigma < 0.0441123197 要求 E > gamma-1+log(2/0.0441123197)。
    E_win = inside(g-1+(2/sig_win_hi).log(), '3.391379', '3.391380')
    E_any = inside(g-1+(2/sig_o).log(), '0.9494525', '0.9494527')
    # 实现：R=2 于 [1,L]、其后为 1。sigma(L)=2e^{gamma-1}(2-1/L)/L 关于 L>1 严格递减，取到 sigma_32。
    sig32 = arb('0.0398046440')                            # 常数锚混合超迹泛函定理 2：sigma_32 in (0.0398046438,0.0398046442)
    for s32, tag in ((arb('0.0398046438'), 'lo'), (arb('0.0398046442'), 'hi')):
        u = 1-(1-s32/floor).sqrt()                          # 2u-u^2 = s32/floor，u=1/L
        Lr = 1/u
        sig, E, _ = step_sigma([arb(1), Lr], [arb(2)])
        assert (sig-s32).contains(0) or abs(sig-s32) < arb('1e-60')
        if tag == 'lo':
            L_hi, E_hi = Lr, E
        else:
            L_lo, E_lo = Lr, E
    inside(L_lo, '65.3397099', '65.3397100')
    inside(L_hi, '65.3397106', '65.3397107')
    inside(E_lo, '3.6834555', '3.6834558')
    inside(E_hi, '3.6834555', '3.6834558')
    inside(L_lo.log(), '4.1795999', '4.1796001')
    full_log = inside(L_hi.log(), '4.1795999', '4.1796001')       # 该例的 int (R-1)_+ ds/s = log L
    assert E_lo > E_win
    print('sigma 下界 2e^(gamma-1)', floor.str(12), '硬截止/热核/最少PV', hard.str(6), heat.str(6), pv.str(7))
    print('进窗所需 E >', E_win.str(7), '有驻点所需 E >=', E_any.str(7), '两阶梯例 L_32', L_lo.str(8), 'E', E_lo.str(5),
          'int(R-1)_+/s', full_log.str(6))
    return floor, L_lo.union(L_hi)


# ---------------------------------------------------------------- 二、四维加六维 CW 混合族
class Sector:
    """T(t) = -K e^{-t} (t-tz)/(t-tp)；Delta=tz-tp。"""

    def __init__(self, K, tz, tp):
        self.K, self.tz, self.tp = K, tz, tp
        self.D = tz-tp
        assert self.D > -4 and not self.D.contains(0)
        if self.D > 0:                                       # 两个临界点 u=(D +- sqrt(D^2+4D))/2，u=t-tp
            r = (self.D*self.D+4*self.D).sqrt()
            self.tm = tp+(self.D-r)/2                         # 在 (-inf,tp) 上，|T| 在此取最小值 m
            self.tpl = tp+(self.D+r)/2                        # 在 (tz,inf) 上，|T| 在此取最大值 M
            assert self.tm < tp and self.tpl > tz
            self.m, self.M = -self.T(self.tm), -self.T(self.tpl)
            assert self.m > 0 and self.M > 0 and self.m > self.M

    def T(self, t):
        return -self.K*(-t).exp()*(t-self.tz)/(t-self.tp)

    def inv(self, tau, a, b):
        """T 在 [a,b] 上严格单调；返回包含 T(t)=tau 之根的球。"""
        fa, fb = self.T(a)-tau, self.T(b)-tau
        assert (fa < 0 < fb) or (fb < 0 < fa), (a, b, tau)
        up = fa < 0
        for _ in range(400):
            if (b-a).upper() < arb(10)**-45:
                break
            mid = arb(((a+b)/2).mid())
            fm = self.T(mid)-tau
            if fm.contains(0):
                break
            if (fm < 0) == up:
                a = mid
            else:
                b = mid
        return a.union(b)

    def comps(self):
        """单调分段（名字, 左端, 右端），端点为搜索用有限数。"""
        eps = arb(10)**-8                                     # tp 的球半径约 1e-10，取远大于它的离极距离
        lo, hi = arb(-80), arb(80)
        if self.D < 0:   # tz<tp：(-inf,tz) 增、(tz,tp) 增、(tp,inf) 增
            return [('C1', lo, self.tz), ('C2', self.tz, self.tp-eps), ('C3', self.tp+eps, hi)]
        return [('D1L', lo, self.tm), ('D1R', self.tm, self.tp-eps), ('D2', self.tp+eps, self.tz),
                ('D3L', self.tz, self.tpl), ('D3R', self.tpl, hi)]

    def pieces(self, ta, tb):
        """负 tau 格 [ta,tb] 内全部根所在的区间（lo, hi），lo/hi 为 None 表示伸向 -inf/+inf（仅尾格出现）。"""
        out = []
        c = {n: (a, b) for n, a, b in self.comps()}
        if self.D < 0:
            for n in ('C1', 'C3'):
                a, b = c[n]
                lo = (None if n == 'C1' else self.tp) if ta is None else self.inv(ta, a, b)
                hi = (self.tz if n == 'C1' else None) if tb is None else self.inv(tb, a, b)
                out.append((n, lo, hi))
            return out
        # D1：谷，T<=-m。根满足 T>=ta，即 |T|<=|ta|。
        if ta is None or ta < -self.m:
            if ta is None:
                L_a, R_a = None, self.tp
            else:
                L_a, R_a = self.inv(ta, *c['D1L']), self.inv(ta, *c['D1R'])
            if tb < -self.m:
                out.append(('D1L', L_a, self.inv(tb, *c['D1L'])))
                out.append(('D1R', self.inv(tb, *c['D1R']), R_a))
            else:
                assert ta is not None
                out.append(('D1', L_a, R_a))
        else:
            assert ta > -self.m
        # D3：峰，-M<=T<0。根满足 T<=tb。
        if tb is None or tb > -self.M:
            if tb is None:
                one, two = self.tz, None
            else:
                one, two = self.inv(tb, *c['D3L']), self.inv(tb, *c['D3R'])
            if ta is not None and ta > -self.M:
                out.append(('D3L', one, self.inv(ta, *c['D3L'])))
                out.append(('D3R', self.inv(ta, *c['D3R']), two))
            else:
                out.append(('D3', one, two))
        else:
            assert tb < -self.M
        return out


def far_ratio(sec, t1):
    """远端 t 在 t1 之外时 r(t)=(t-tz)/(t-tp)=1+D/(t-tp) 介于 1 与 r(t1) 之间。"""
    r1 = 1+sec.D/(t1-sec.tp)
    assert r1 > 0
    return (r1.log()).union(arb(0))


def diff_range(so, se, po, pe):
    """(t_e - t_o) 的包围；两端同向伸出时用远端渐近 t=log(K r/|tau|)。"""
    _, olo, ohi = po
    _, elo, ehi = pe
    if olo is None and elo is None:                           # 同向趋于 -inf
        return (se.K/so.K).log()+far_ratio(se, ehi)-far_ratio(so, ohi), 'far-'
    if ohi is None and ehi is None:                           # 同向趋于 +inf
        return (se.K/so.K).log()+far_ratio(se, elo)-far_ratio(so, olo), 'far+'
    lo = None if (elo is None or ohi is None) else elo-ohi
    hi = None if (ehi is None or olo is None) else ehi-olo
    return (lo, hi), 'box'


def excluded(rng, kind, Wlo, Whi):
    if kind != 'box':
        return rng < Wlo or rng > Whi
    lo, hi = rng
    return (hi is not None and hi < Wlo) or (lo is not None and lo > Whi)


def cw46(S, L, Wlo, Whi):
    results = {}
    for name, delta in (('0', arb(0)), ('1/3', arb(1)/3)):
        sec = {}
        for s in 'oe':
            K = 2*S[(s, 2)]/(3*S[(s, 3)])
            assert K > 0
            tz = -L[(s, 2)]/S[(s, 2)]
            tp = -L[(s, 3)]/S[(s, 3)]+delta+arb(1)/6
            sec[s] = Sector(K, tz, tp)
        so, se = sec['o'], sec['e']
        assert se.D < 0                                       # 两个 delta 下偶荷都是三段单增
        assert (so.D > 0) == (name == '0')
        # tau>0：奇、偶各恰一根，分别在 tz 与 tp 之间，k 有界。
        lo_o, hi_o = so.tz.min(so.tp), so.tz.max(so.tp)
        lo_e, hi_e = se.tz.min(se.tp), se.tz.max(se.tp)
        kpos = ((lo_e-hi_o)/2).exp(), ((hi_e-lo_o)/2).exp()
        assert (hi_e-lo_o) < Wlo
        # tau<0：自适应格。中段格的端点都在 [-10^4,-10^-5] 内；两端尾格各一。
        cand = []
        tmin, tmax = arb(10)**-5, arb(10)**4
        # 尾格 (-inf,-tmax]
        po_all = so.pieces(None, -tmax)
        pe_all = se.pieces(None, -tmax)
        for po in po_all:
            for pe in pe_all:
                rng, kind = diff_range(so, se, po, pe)
                assert excluded(rng, kind, Wlo, Whi), ('tail-', name, po, pe, rng)
        # 尾格 [-tmin,0)：根在 tz 附近或伸向 +inf
        po_all = so.pieces(-tmin, None)
        pe_all = se.pieces(-tmin, None)
        for po in po_all:
            for pe in pe_all:
                rng, kind = diff_range(so, se, po, pe)
                assert excluded(rng, kind, Wlo, Whi), ('tail+', name, po, pe, rng)
        # 中段：log|tau| 从 4log10 到 -5log10，先均分，再对不能排除的格二分。
        N = 180
        lt = [arb(4)*arb(10).log() - arb(9)*arb(10).log()*arb(i)/N for i in range(N+1)]
        stack = [(-(lt[i].exp()), -(lt[i+1].exp()), 0) for i in range(N)]
        stack = [(arb(a.mid()), arb(b.mid()), dpt) for a, b, dpt in stack]
        stack[0] = (-tmax, stack[0][1], 0)
        stack[-1] = (stack[-1][0], -tmin, 0)
        while stack:
            ta, tb, dpt = stack.pop()
            assert ta < tb
            ok, interior = True, False
            for po in so.pieces(ta, tb):
                for pe in se.pieces(ta, tb):
                    rng, kind = diff_range(so, se, po, pe)
                    if not excluded(rng, kind, Wlo, Whi):
                        ok = False
                        if kind == 'box' and None not in rng and rng[0] > Wlo and rng[1] < Whi:
                            interior = True               # 整格都落在窗内的格不必再细分
            if ok:
                continue
            if interior or dpt >= 18:
                cand.append((ta, tb))
                continue
            mid = arb(((ta+tb)/2).mid())
            stack += [(ta, mid, dpt+1), (mid, tb, dpt+1)]
        cand.sort(key=lambda c: c[0].mid())
        merged = []
        for a, b in cand:
            if merged and (merged[-1][1]-a).contains(0) or (merged and merged[-1][1] >= a):
                merged[-1] = (merged[-1][0], b)
            else:
                merged.append((a, b))
        results[name] = (sec, merged, kpos)
    return results


def natural_members(sec, taus, Whi):
    """固定 tau 的成员：逐段求出全部根，列出全部驻点倍率。"""
    so, se = sec['o'], sec['e']
    out = {}
    for tau in taus:
        roots = {}
        for s, sc in (('o', so), ('e', se)):
            rs = []
            for n, a, b in sc.comps():
                fa, fb = sc.T(a)-tau, sc.T(b)-tau
                if (fa < 0 < fb) or (fb < 0 < fa):
                    rs.append(sc.inv(tau, a, b))
                else:
                    assert (fa > 0 and fb > 0) or (fa < 0 and fb < 0)   # 单调段两端同号：无根
            roots[s] = rs
        ks = [((te-to)/2).exp() for to in roots['o'] for te in roots['e']]
        top = None
        for k in ks:
            top = k if top is None else top.max(k)
        out[tau.str(5)] = (ks, top)
    return out


def verify(root):
    S, L, hmax2 = sector_sums(root)
    A = {s: S[(s, 1)]/S[(s, 2)]*(L[(s, 2)]/S[(s, 2)]).exp() for s in 'oe'}
    sig = {s: 2/(arb(1).exp()*A[s]) for s in 'oe'}
    inside(sig['o'], '0.50707835', '0.50707836')
    inside(sig['e'], '0.6683849', '0.6683851')
    floor, L32 = schwinger_part(sig['o'], sig['e'], arb('0.0441123197'))
    # 窗内截断驻点处的 x_r/M0 = 2 h_r^2 e^{w_s - L4_s/S4_s}/sigma 只依赖 sigma（对 S_+ 中任意 R）。
    s32 = arb('0.0398046438').union(arb('0.0398046442'))
    off = {s: L[(s, 2)]/S[(s, 2)] for s in 'oe'}
    wo = (-s32*A['o']/2).lambertw(branch=-1)
    we = (-s32*A['e']/2).lambertw(branch=0)
    assert wo.is_finite() and we.is_finite()
    ratio_e = inside(2*hmax2['e']*(we-off['e']).exp()/s32, '94.8', '95.0')
    ratio_o = inside(2*hmax2['o']*(wo-off['o']).exp()/s32, '0.36', '0.38')
    # 两阶梯例 R=2 于 [1,L]：M0=2-1/L，C=1-gamma+log L。
    M0 = inside(2-1/L32, '1.98', '1.99')
    C = 1-arb.const_euler()+L32.log()
    xe = inside(ratio_e*M0, '188.1', '188.5')
    xo = inside(ratio_o*M0, '0.72', '0.75')
    fT = inside(M0+xe*xe.log()-C*xe, '120', '122')
    assert fT > 60*M0                                        # 精确势 0<f'<M0≈1.98，截断在那里不逼近精确势
    print('sigma_32 处最重偶/奇荷 x/M0', ratio_e.str(5), ratio_o.str(3), '两阶梯例 x_max 偶/奇', xe.str(5), xo.str(3),
          '截断 f_T\'(x_e)', fT.str(4), 'M0', M0.str(6))

    klo, khi = arb('29.9841147'), arb('35.0010644')
    Wlo, Whi = 2*klo.log(), 2*khi.log()
    k0 = ((L[('o', 2)]/S[('o', 2)]-L[('e', 2)]/S[('e', 2)])/2).exp()
    inside(k0, '2.4134605455', '2.4134605465')
    res = cw46(S, L, Wlo, Whi)
    expect = {'0': [('-0.0042914', '-0.0031447'), ('-108.066', '-79.338')],
              '1/3': [('-0.0046164', '-0.0033659'), ('-72.345', '-52.846')]}
    for name, (sec, merged, kpos) in res.items():
        so, se = sec['o'], sec['e']
        print('delta', name, 'Delta_o', so.D.str(6), 'Delta_e', se.D.str(6),
              'tau>0 时 k 在', kpos[0].str(6), kpos[1].str(6))
        assert len(merged) == 2, merged
        ivs = sorted(merged, key=lambda c: c[0].mid())
        big, small = ivs[0], ivs[1]
        print('  进窗只可能在 tau 属于', small[0].str(7), small[1].str(7), '或', big[0].str(6), big[1].str(6))
        (s_lo, s_hi), (b_lo, b_hi) = expect[name]
        assert arb(s_lo) < small[0] and small[1] < arb(s_hi)
        assert arb(b_lo) < big[0] and big[1] < arb(b_hi)
        # 每个候选格并集内确有进窗点：取中点逐段求根。
        for a, b in (small, big):
            m = (a+b)/2
            mem = natural_members(sec, [m], Whi)
            ks = list(mem.values())[0][0]
            assert any(klo < k < khi for k in ks), (name, m)
        if name == '0':
            assert so.m > 0 and so.M > 0
            inside(so.m, '0.9519', '0.9520')
            inside(so.M, '0.37266', '0.37267')
        mem = natural_members(sec, [arb(1).exp(), -arb(1).exp(), arb(1), -arb(1)], Whi)
        for t, (ks, top) in mem.items():
            print('  tau', t, '驻点数', len(ks), '最大倍率', top.str(6) if top is not None else None)
            assert top is not None and top < klo
        want = {'0': ('5.5673', '5.5674', '3.3516', '3.3517'), '1/3': ('7.3809', '7.3810', '5.1617', '5.1618')}[name]
        inside(mem[(-arb(1).exp()).str(5)][1], want[0], want[1])
        inside(mem[(-arb(1)).str(5)][1], want[2], want[3])
        pos = {'0': (('2.60416', '2.60418'), ('2.55742', '2.55744')), '1/3': (('2.61627', '2.61629'), ('2.57157', '2.57159'))}[name]
        inside(mem[arb(1).exp().str(5)][1], *pos[0])
        inside(mem[arb(1).str(5)][1], *pos[1])
        assert [len(mem[t.str(5)][0]) for t in (arb(1).exp(), -arb(1).exp(), arb(1), -arb(1))] == [1, 4, 1, 4]
        kp = {'0': ('2.41345', '2.64734'), '1/3': ('2.09961', '3.04305')}[name]
        assert arb(kp[0]) < kpos[0] and kpos[1] < arb(kp[1])
        # 原像相对宽度（端点比）
        print('  两区间端点比', (small[0]/small[1]).str(6), (big[0]/big[1]).str(6), 'lambda=tau/e（c4=3/2）',
              (small[0]/arb(1).exp()).str(5), (small[1]/arb(1).exp()).str(5), (big[0]/arb(1).exp()).str(5), (big[1]/arb(1).exp()).str(5))
        rat = {'0': ('1.36457', '1.36459', '1.36206', '1.36208'), '1/3': ('1.37147', '1.37149', '1.36894', '1.36896')}[name]
        inside(small[0]/small[1], rat[0], rat[1])
        dd = {'0': (('0.0547250', '0.0547252'), ('-0.130258', '-0.130256')),
              '1/3': (('-0.278609', '-0.278607'), ('-0.463591', '-0.463589'))}[name]
        inside(so.D, *dd[0])
        inside(se.D, *dd[1])
        lam = {'0': ('-0.001579', '-0.001156', '-39.76', '-29.18'), '1/3': ('-0.0016983', '-0.0012382', '-26.615', '-19.44')}[name]
        e1 = arb(1).exp()
        assert arb(lam[0]) < small[0]/e1 and small[1]/e1 < arb(lam[1])
        assert arb(lam[2]) < big[0]/e1 and big[1]/e1 < arb(lam[3])
        inside(big[0]/big[1], rat[2], rat[3])
        # 独立核对（浮点直接扫描驻定方程 G_s(t)=0 的变号数，不经 T_s 的分段）：
        import math
        for tau in (math.e, -math.e, 1.0, -1.0, float(((small[0]+small[1])/2).mid()), float(((big[0]+big[1])/2).mid())):
            cnt = {}
            for s_ in 'oe':
                S4, S6 = float(S[(s_, 2)].mid()), float(S[(s_, 3)].mid())
                L4, L6 = float(L[(s_, 2)].mid()), float(L[(s_, 3)].mid())
                dl = 0.0 if name == '0' else 1/3
                def G(t):
                    return 2*(S4*t+L4)+tau*math.exp(t)*(3*S6*t+3*L6-(3*dl+0.5)*S6)
                grid = [-60+i*0.0005 for i in range(160001)]
                vals = [G(t) for t in grid]
                cnt[s_] = sum(1 for a, b in zip(vals, vals[1:]) if a*b < 0)
            ref = natural_members(sec, [arb(tau)], Whi)
            assert cnt['o']*cnt['e'] == len(list(ref.values())[0][0]), (name, tau, cnt)
    return res


if __name__ == '__main__':
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
