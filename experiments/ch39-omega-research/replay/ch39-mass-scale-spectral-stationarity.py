"""两个质量标尺 (k,d) 上的单谱型泛函驻点与沿 k/d^2 滑行的读数响应。

记 u=log k, v=log d。奇荷 (1,3) 重模质量为 Mbar*h/d^2，偶荷 (0,2,4,6) 为 Mbar*k*h/d^2，
故全部质量只经两条线性型 alpha=-2v 与 beta=u-2v 依赖标尺；沿 k/d^2 滑行即 beta 固定、alpha 变动。
"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


FLAGSHIP = 2_306_888_100_085


def field_rows(vacuum, layers):
    """与原 minimal['fields'] 同序的 (g, 牛顿权, 偶荷?, h)，并逐行核对。"""
    rows = []
    for q, dim in ((1, 6), (4, 3), (2, 3), (3, 2), (6, 1)):
        for h in layers[q][3:]:
            extra = 4 if q == 3 else 0
            rows.append((-4*dim+extra, arb(dim)/3+arb(extra)/6, q % 2 == 0, h))
    rows.extend((36, arb(-6), True, h) for h in layers[0])
    original = vacuum['minimal']['fields'](layers)
    assert len(original) == len(rows)+1 and original[0][0] == -62 and original[0][2].is_zero()
    scale = arb('2.435e18')/(724**2*arb('1.22089e19'))
    for (g, a, z), (g2, a2, even, h) in zip(original[1:], rows, strict=True):
        assert g == g2 and (a-a2).contains(0)
        assert (z-((32 if even else 1)*scale*h)**2).contains(0)
    return rows


def weighted_log_mean(pairs):
    den = sum(w for w, _ in pairs)
    assert not den.contains(0)
    return sum(w*h.log() for w, h in pairs)/den, den


def cw_point(rows):
    """四维一圈 STr M^4(log M^2/Lambda^2 - c) 的驻点倍率 k* 与奇荷对数均值。"""
    odd, ode = weighted_log_mean([(g*h**4, h) for g, _, e, h in rows if not e])
    even, eve = weighted_log_mean([(g*h**4, h) for g, _, e, h in rows if e])
    return (odd-even).exp(), odd, ode, eve


def layers_at(chiral, model, n):
    if n < 10000:
        # 与 ch39-first-three-mode-observables.py 相同的按编号重建。
        levels = model['shell']['shell_spectra']((0, *model['charges']), zero_index=n)
        center = sum(levels[0])/18
        return {q: chiral['ordered']([v*abs(v-center) for v in eig]) for q, eig in levels.items()}
    rows = chiral['source']['spectral_supplier'](model)(acb.zeta_zero(n).imag)
    return {q: chiral['ordered']([h for h, z in rows if z == q]) for q in (0, *model['charges'])}


def verify(root):
    physical = run_path(str(root/'ch39-first-three-mode-observables.py'))
    model, initial = physical['model_inputs']()
    ctx.prec, ctx.cap = 192, 3
    joint = physical['joint']
    chiral, vacuum = joint['chiral'], joint['vacuum']
    layers = layers_at(chiral, model, FLAGSHIP)
    assert all(layers[q][2] < 9 < layers[q][3] for q in model['charges'])
    rows = field_rows(vacuum, layers)

    # 一、整数计数：奇荷行全为费米型负权，偶荷行混号。
    odd_g = [g for g, _, e, _ in rows if not e]
    even_g = [g for g, _, e, _ in rows if e]
    assert len(odd_g) == 32 and len(even_g) == 66
    assert set(odd_g) == {-24, -4} and sum(odd_g) == -448 and sum(even_g) == 200
    assert sum(odd_g)+sum(even_g)-62 == -310
    # STr log M^2 对 (u,v) 的梯度为常向量 (2G_e, -4(G_o+G_e))。
    gradient = (2*sum(even_g), -4*(sum(odd_g)+sum(even_g)))
    assert gradient == (400, 992)
    assert 2*gradient[0]+gradient[1] == 1792  # 沿滑行方向 (du,dv)=(2,1) 的导数
    # 四层低能区（degree=4）按 low_sector 的定义取 g=s-4f；其质量平方 j*a、a ∝ (k/d^2)^-2。
    fermions = (3*4+3)//4
    low_g = [s-4*f for s, f in zip([-4, 12, 4*fermions-12, 4], [0, 0, fermions, 0])]
    assert low_g == [-4, 12, -12, 4] and sum(low_g[1:]) == 4
    rows_low = vacuum['low_sector'](4, arb(1))
    assert all((z-j).is_zero() for j, (_, z) in enumerate(rows_low))
    # d log(j a)/d(u,v) = (-2, 4)，故梯度加上 4*(-2,4)。
    full_gradient = (gradient[0]+4*(-2), gradient[1]+4*4)
    assert full_gradient == (392, 1008) and 2*full_gradient[0]+full_gradient[1] == 1792
    # 规范交点只读奇荷层：改动全部偶荷层，阈值表逐项不变。
    two = joint['gauge']['two']
    shifted = {q: ([2*h for h in band] if q % 2 == 0 else band) for q, band in layers.items()}
    table, table2 = (two['thresholds'](model, x, False) for x in (layers, shifted))
    assert len(table) == len(table2) == 32
    assert all(a[0].mid() == b[0].mid() and a[3:] == b[3:] and str(a[1:3]) == str(b[1:3])
               for a, b in zip(table, table2))
    # 牛顿热核权及原配平后的残余权在奇荷行上同样严格同号。
    for g, a, e, _ in rows:
        if not e:
            residual = a+(arb(-g)/6 if g < 0 else arb(g)/12)
            assert a > 0 and residual > 0
    print('奇/偶荷行数与超迹计数', len(odd_g), len(even_g), sum(odd_g), sum(even_g),
          'STr log 梯度', gradient)

    # 二、四维 Coleman–Weinberg 型：驻点倍率与锚标度 Lambda、常数 c 无关。
    kstar, odd_mean, odd_den, even_den = cw_point(rows)
    assert odd_den < 0 and even_den > 0
    assert arb('2.4134') < kstar < arb('2.4135'), kstar
    # 驻点 Hessian：奇向 8 Sum_o g m^4 < 0，偶向 8 Sum_e g m^4 > 0，交叉项为零，故为鞍点。
    # 奇荷驻点：log d*^2 = log(Mbar/Lambda) + <log h>_o - (c-1/2)/2。
    dstar = (odd_mean-arb(1)/2).exp().sqrt()
    assert arb(9) < dstar < arb('9.6')
    lam_724 = arb('2.435e18')*(odd_mean-arb(1)/2-2*arb(724).log()).exp()
    print('CW 驻点倍率 k*', kstar.str(12), '奇荷 h^4|g| 对数均值', odd_mean.str(12))
    print('Lambda=Mbar, c=3/2 的 d*', dstar.str(10), '使 d*=724 所需 Lambda/GeV', lam_724.str(8))

    # 三、任意维数 D=2p 的 CW 型，倍数权与 |g| 权：对全部 p>0 的上下界。
    hmax_o = max(h for _, _, e, h in rows if not e)
    hmax_e = max(h for _, _, e, h in rows if e)
    bounds = {}
    for name, weight in (('倍数', lambda g: 1), ('|g|', abs)):
        geo_o, _ = weighted_log_mean([(arb(weight(g)), h) for g, _, e, h in rows if not e])
        geo_e, _ = weighted_log_mean([(arb(weight(g)), h) for g, _, e, h in rows if e])
        bounds[name] = (geo_o.exp()/hmax_e, hmax_o/geo_e.exp())
        print('全部 p>0 的 k*(p) 包含于', name, *(v.str(8) for v in bounds[name]))

    # 四、已认证两窗与 k>=1 给出的实倍率外区间 (证书见 ch39-joint-integer-scale-pairs.py)。
    klo = arb('31.98396690')*(arb(701)/724)**2
    khi = arb('32.01596692')*(arb(757)/724)**2
    assert arb('29.98') < klo < arb('29.99') and arb('35.00') < khi < arb('35.01')
    assert kstar < klo and all(b[1] < klo for b in bounds.values())
    print('两窗允许的实倍率外区间', klo.str(8), khi.str(8))
    # 逐行常数 c_r 落在长为 Delta 的区间时 |log k* - log k*_0| <= Delta (R_e+1)/4。
    ratio = sum(abs(g)*h**4 for g, _, e, h in rows if e)/even_den
    delta_max = 4*(klo/kstar).log()/(ratio+1)
    assert arb('10.2') < ratio < arb('10.3') and delta_max > arb('0.89') > arb(2)/3
    print('偶荷 |g|h^4 与 g h^4 之比', ratio.str(8), '容许的逐行常数跨度上界', delta_max.str(8))

    # 五、另两枚实际零点上的同一驻点倍率。
    others = {}
    for n, lo, hi in ((11661, '2.3', '2.6'), (8359, '2.3', '2.6')):
        rows_n = field_rows(vacuum, layers_at(chiral, model, n))
        k_n, _, od, ev = cw_point(rows_n)
        assert od < 0 and ev > 0 and arb(lo) < k_n < arb(hi)
        others[n] = k_n
        print('零点', n, 'CW 驻点倍率', k_n.str(10))

    # 六、沿 k/d^2=32/724^2 的滑行：偶质量、跷跷板与曲率主项不变，只有奇阈值移动。
    k0, d0 = arb(32), arb(724)
    mp, reduced, vev = map(arb, ('1.22089e19', '2.435e18', '246.22'))
    coefficient = (9*arb(3).log()-12*arb(2).log())/2
    a0 = (vev**2/(mp*(k0*reduced/d0**2)*layers[2][15]))**2
    main = coefficient*a0**2/(arb.pi()*arb('1e-122'))
    loss = 2*a0/(arb('.006')*coefficient)
    assert loss < arb('1e-57') and arb('.998') < main*(1-loss) and main < arb('1.002')
    even_masses = [*layers[0], *(h for q in (2, 4, 6) for h in layers[q][3:])]
    lightest_even = k0*reduced/d0**2*min(even_masses)
    strong = {}
    for d in (715, 716, 724, 743, 744):
        dd = arb(d)
        bands = {q: [h*(d0/dd)**2*((k0*(dd/d0)**2)/k0 if q % 2 == 0 else 1) for h in band]
                 for q, band in layers.items()}
        # 偶层精确不变：k=32(d/724)^2 时 (d0/d)^2 (k/k0) = 1。
        assert all((x-y).contains(0) for q in (0, 2, 4, 6) for x, y in zip(bands[q], layers[q]))
        residual = [(w+(arb(-g)/6 if g < 0 else arb(g)/12), z)
                    for g, w, z in vacuum['minimal']['fields'](bands)]
        residual += vacuum['low_sector'](4, a0)
        assert vacuum['newton'](residual, arb('.03')) < 1 < vacuum['newton'](residual, arb('.033'))
        with redirect_stdout(StringIO()):
            energy, _, alpha = joint['gauge']['crossing_enclosure'](
                model, initial, layers, False, dd, '29.04', '.1184', recenter=True)
        assert lightest_even > energy
        strong[d] = (energy, alpha)
        print('滑行点 d =', d, 'k =', (k0*(dd/d0)**2).str(12), '交点/GeV', energy.str(8),
              '强耦合', alpha.str(10))
    window = lambda x: arb('.117') < x < arb('.119')
    assert window(strong[716][1]) and window(strong[724][1]) and window(strong[743][1])
    assert strong[715][1] > arb('.119') and strong[744][1] < arb('.117')
    secant = (strong[743][1]-strong[716][1])/27
    assert arb('-7.05e-5') < secant < arb('-6.95e-5')
    log_secant = (strong[743][1]-strong[716][1])/(arb(743).log()-arb(716).log())
    assert arb('-0.0515') < log_secant < arb('-0.0505')
    print('716 至 743 的强耦合割线 /单位分母、/单位 log d', secant.str(8), log_secant.str(8))
    print('曲率主项（滑行不变）', main.str(14), '修正上界', loss.str(4))


if __name__ == '__main__':
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
