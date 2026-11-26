"""两个质量标尺 (k,d) 上的超迹谱泛函驻点（整数维与实数维）与沿 k/d^2 滑行的读数响应。

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


def inside(x, lo, hi):
    assert arb(lo) < x < arb(hi), (x, lo, hi)
    return x


def weighted_log_mean(pairs):
    den = sum(w for w, _ in pairs)
    assert not den.contains(0)
    return sum(w*h.log() for w, h in pairs)/den, den


def power_point(odd, even, p):
    """D=2p 维 STr M^{2p}(log M^2/Lambda^2-c) 型的驻点倍率；odd/even 为 (权, h)。"""
    so = sum(g*h**(2*p) for g, h in odd)
    se = sum(g*h**(2*p) for g, h in even)
    lo = sum(g*h**(2*p)*h.log() for g, h in odd)/so
    le = sum(g*h**(2*p)*h.log() for g, h in even)/se
    return (lo-le).exp(), lo, so, se


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
    # 四层低能区（degree=4）：独立写出的场变化字面量，g=s-4f；与原模块返回的牛顿权与质量平方逐项比较。
    scalars, diracs = [-4, 12, 0, 4], [0, 0, 3, 0]
    low_g = [s_-4*f for s_, f in zip(scalars, diracs)]
    assert low_g == [-4, 12, -12, 4] and sum(low_g[1:]) == 4
    rows_low = vacuum['low_sector'](4, arb(1))
    assert len(rows_low) == 4
    for j, ((w, z), s_, f) in enumerate(zip(rows_low, scalars, diracs)):
        assert (z-j).is_zero() and (w-(arb(s_)/6+arb(f)/3)).contains(0)
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

    # 二、四维 Coleman–Weinberg 超迹：驻点倍率与锚标度 Lambda、常数 c 无关。
    odd_w = [(arb(g), h) for g, _, e, h in rows if not e]
    even_w = [(arb(g), h) for g, _, e, h in rows if e]
    kstar, odd_mean, odd_den, even_den = power_point(odd_w, even_w, 2)
    assert odd_den < 0 < even_den
    inside(kstar, '2.4134605455', '2.4134605465')
    inside(odd_mean, '5.00398096', '5.00398097')
    # Hessian：取 Lambda=Mbar、c=3/2，以 e^alpha 为奇荷质量单位直接验算驻点与二阶导数。
    # W(alpha)=Sum g m^4 (log m^2/Lambda^2 - c)，W'=Sum g m^4 [4(l-c)+2]，W''=Sum g m^4 [16(l-c)+16]。
    c0 = arb(3)/2
    for sector, mean in ((odd_w, odd_mean), (even_w, None)):
        if mean is None:
            mean = sum(g*h**4*h.log() for g, h in sector)/sum(g*h**4 for g, h in sector)
        scale_log = (c0-arb(1)/2)/2-mean          # 驻点处 alpha-log(Lambda)
        first = sum(g*((scale_log+h.log())*4).exp()*(4*(2*(scale_log+h.log())-c0)+2) for g, h in sector)
        second = sum(g*((scale_log+h.log())*4).exp()*(16*(2*(scale_log+h.log())-c0)+16) for g, h in sector)
        size = sum(abs(g)*((scale_log+h.log())*4).exp() for g, h in sector)
        assert first.contains(0) and abs(first) < size*arb('1e-8')
        assert (second < 0) if sector is odd_w else (second > 0)  # 奇向极大、偶向极小：鞍点
    # 奇荷驻点：log d*^2 = log(Mbar/Lambda) + <log h>_o - (c-1/2)/2。
    dstar = inside((odd_mean-arb(1)/2).exp().sqrt(), '9.506639', '9.506641')
    lam_724 = inside(arb('2.435e18')*(odd_mean-arb(1)/2-2*arb(724).log()).exp(), '4.198323e14', '4.198324e14')
    print('CW 驻点倍率 k*', kstar.str(12), '奇荷 g h^4 对数均值', odd_mean.str(12))
    print('Lambda=Mbar, c=3/2 的 d*', dstar.str(10), '使 d*=724 所需 Lambda/GeV', lam_724.str(8))

    # 三、两窗与 k>=1 给出的实倍率外区间。曲率必要带按 ch39-joint-integer-scale-pairs.py 的公式重算，
    # 并核对该证书所用的外端字面量确实包住它；其实平移排除覆盖 s>=-3.5 与 d 不在 (701,757)。
    source = (root/'ch39-joint-integer-scale-pairs.py').read_text(encoding='utf-8')
    kmin_used, kmax_used = arb('31.98396690'), arb('32.01596692')
    assert "arb('31.98396690'), arb('32.01596692')" in source
    k0, d0 = arb(32), arb(724)
    mp, reduced, vev = map(arb, ('1.22089e19', '2.435e18', '246.22'))
    coefficient = (9*arb(3).log()-12*arb(2).log())/2
    a0 = (vev**2/(mp*(k0*reduced/d0**2)*layers[2][15]))**2
    main = coefficient*a0**2/(arb.pi()*arb('1e-122'))
    band_lo = 32*(main*(1-arb('1e-57'))/arb('1.002')).root(4)
    band_hi = 32*(main/arb('.998')).root(4)
    assert kmin_used < band_lo and band_hi < kmax_used
    klo = inside(kmin_used*(arb(701)/724)**2, '29.9841147', '29.9841153')
    khi = inside(kmax_used*(arb(757)/724)**2, '35.0010638', '35.0010644')
    dmin = inside(724/kmax_used.sqrt(), '127.95', '127.96')  # k>=1 时 d 的下界
    assert 2*(dmin/724).log() > arb('-3.5')
    print('曲率必要带', band_lo.str(12), band_hi.str(12), 'k>=1 时 d >', dmin.str(8))
    print('两窗允许的实倍率外区间（k>=1）', klo.str(10), khi.str(10))
    assert kstar < klo
    # 逐行常数 c_r 落在长为 Delta 的区间时 |log k* - log k*_0| <= Delta (R_e+1)/4。
    ratio = inside(sum(abs(g)*h**4 for g, h in even_w)/even_den, '10.212445', '10.212447')
    delta_max = inside(4*(klo/kstar).log()/(ratio+1), '0.8988604', '0.8988605')
    print('偶荷 |g|h^4 与 g h^4 之比', ratio.str(8), '容许的逐行常数跨度上界', delta_max.str(8))

    # 四、D=2p 的有符号超迹族。
    # (a) 全部整数 D=1..59 逐个计算；(b) p>=30 的实数尾部由最高偶模主导给出解析包络。
    values = {}
    for D in range(1, 60):
        k_d, _, so, se = power_point(odd_w, even_w, arb(D)/2)
        assert so < 0 and (se > 0 if D <= 11 else se < 0)
        assert k_d < klo
        values[D] = k_d
    inside(values[11], '15.5470', '15.5472')
    inside(values[12], '0.283799', '0.283801')
    assert all(values[D] < values[11] for D in values if D != 11)
    assert all(values[D] > values[12] for D in values if D != 12)
    top = max(even_w, key=lambda t: t[1].mid())
    rest = [t for t in even_w if t is not top]
    assert top[0] == -4 and all(h < top[1]*arb('0.95') for _, h in rest)
    P, H = 30, top[1]
    eps = sum(abs(g)*(h/H)**(2*P) for g, h in rest)
    eta = sum(abs(g)*(H/h).log()*(h/H)**(2*P) for g, h in rest)
    assert eps < 1
    tau = eta/(4-eps)
    _, odd_floor, _, _ = power_point(odd_w, even_w, P)
    hmax_o = max(h for _, h in odd_w)
    tail_lo = inside((odd_floor-H.log()-tau).exp(), '1.9562', '1.9564')
    tail_hi = inside(hmax_o/H*tau.exp(), '2.0274', '2.0276')
    assert tail_hi < klo
    print('整数 D=1..59 的 k* 最大值（D=11）与最小值（D=12）', values[11].str(8), values[12].str(8))
    print('p>=30 的全部实数 k*(p) 包含于', tail_lo.str(8), tail_hi.str(8))
    # (c) 实数 D 上可取到 32：偶荷分母在 (5.75,5.76) 变号，而 k*(5.55)<32<k*(5.58) 且其间分母为正。
    assert power_point(odd_w, even_w, arb('5.75'))[3] > 0 > power_point(odd_w, even_w, arb('5.76'))[3]
    pieces = 60
    for i in range(pieces):
        ball = arb('5.55')+arb('.03')*i/pieces
        ball = ball.union(arb('5.55')+arb('.03')*(i+1)/pieces)
        assert sum(g*h**(2*ball) for g, h in even_w) > 0
        assert sum(g*h**(2*ball) for g, h in odd_w) < 0
    k_left = power_point(odd_w, even_w, arb('5.55'))[0]
    k_right = power_point(odd_w, even_w, arb('5.58'))[0]
    assert k_left < 32 < k_right
    print('实维数 D=2p：k*(5.55), k*(5.58)', k_left.str(8), k_right.str(8))

    # 五、非超迹的正权族（行倍数权与 |g| 权）：对全部 p>0 的包络。
    hmax_e = max(h for _, h in even_w)
    expected = {'倍数': (('0.4862600', '0.4862601'), ('3.940431', '3.940432')),
                '|g|': (('0.4841085', '0.4841086'), ('4.323695', '4.323696'))}
    for name, weight in (('倍数', lambda g: 1), ('|g|', abs)):
        geo_o, _ = weighted_log_mean([(arb(weight(g)), h) for g, _, e, h in rows if not e])
        geo_e, _ = weighted_log_mean([(arb(weight(g)), h) for g, _, e, h in rows if e])
        low_end = inside(geo_o.exp()/hmax_e, *expected[name][0])
        high_end = inside(hmax_o/geo_e.exp(), *expected[name][1])
        assert high_end < klo
        print('正权族全部 p>0 的 k*(p) 包含于', name, low_end.str(8), high_end.str(8))

    # 六、另两枚实际零点上的四维 CW 驻点倍率。
    for n, lo, hi in ((11661, '2.3788110865', '2.3788110875'), (8359, '2.3830181665', '2.3830181675')):
        rows_n = field_rows(vacuum, layers_at(chiral, model, n))
        k_n, _, od, ev = power_point([(arb(g), h) for g, _, e, h in rows_n if not e],
                                     [(arb(g), h) for g, _, e, h in rows_n if e], 2)
        assert od < 0 < ev
        inside(k_n, lo, hi)
        print('零点', n, 'CW 驻点倍率', k_n.str(10))

    # 七、沿 k/d^2=32/724^2 的滑行：逐点重算偶层、跷跷板 a、曲率主项与牛顿根；只有奇阈值移动。
    inside(main, '0.99999336383', '0.99999336387')
    even_masses = [*layers[0], *(h for q in (2, 4, 6) for h in layers[q][3:])]
    strong = {}
    expected_energy = {715: ('3.85045e14', '3.85109e14'), 716: ('3.83667e14', '3.83755e14'),
                       724: ('3.72980e14', '3.73062e14'), 743: ('3.49022e14', '3.49102e14'),
                       744: ('3.47811e14', '3.47890e14')}
    expected_alpha = {715: ('0.1190155', '0.1190203'), 716: ('0.1189429', '0.1189493'),
                      724: ('0.1183751', '0.1183814'), 743: ('0.1170532', '0.1170593'),
                      744: ('0.1169842', '0.1169904')}
    for d in (715, 716, 724, 743, 744):
        dd = arb(d)
        kk = k0*(dd/d0)**2
        bands = {q: [h*(d0/dd)**2*(kk/k0 if q % 2 == 0 else 1) for h in band]
                 for q, band in layers.items()}
        assert all((x-y).contains(0) for q in (0, 2, 4, 6) for x, y in zip(bands[q], layers[q]))
        a = (vev**2/(mp*(kk*reduced/dd**2)*bands[2][15]))**2
        assert (a-a0).contains(0)
        curve_main = coefficient*a**2/(arb.pi()*arb('1e-122'))
        loss = 2*a/(arb('.006')*coefficient)
        assert loss < arb('8.50e-59')
        assert arb('.998') < curve_main*(1-loss) and curve_main < arb('1.002')
        residual = [(w+(arb(-g)/6 if g < 0 else arb(g)/12), z)
                    for g, w, z in vacuum['minimal']['fields'](bands)]
        residual += vacuum['low_sector'](4, a)
        assert vacuum['newton'](residual, arb('.03')) < 1 < vacuum['newton'](residual, arb('.033'))
        with redirect_stdout(StringIO()):
            energy, _, alpha = joint['gauge']['crossing_enclosure'](
                model, initial, layers, False, dd, '29.04', '.1184', recenter=True)
        assert kk*reduced/dd**2*min(even_masses) > energy
        inside(energy, *expected_energy[d])
        inside(alpha, *expected_alpha[d])
        assert alpha.rad() < arb('3.3e-6')
        strong[d] = alpha
        print('滑行点 d =', d, 'k =', kk.str(12), '交点/GeV', energy.lower().str(7), energy.upper().str(7),
              '强耦合', alpha.lower().str(8), alpha.upper().str(8), '曲率修正上界', loss.str(4))
    window = lambda x: arb('.117') < x < arb('.119')
    assert window(strong[716]) and window(strong[724]) and window(strong[743])
    assert strong[715] > arb('.119') and strong[744] < arb('.117')
    secant = inside((strong[743]-strong[716])/27, '-7.03e-5', '-6.97e-5')
    assert secant.rad() < arb('2.4e-7')
    log_secant = inside((strong[743]-strong[716])/(arb(743).log()-arb(716).log()), '-0.0513', '-0.0507')
    print('716 至 743 的强耦合割线 /单位分母、/单位 log d', secant.lower().str(6), secant.upper().str(6),
          log_secant.lower().str(6), log_secant.upper().str(6))
    print('曲率主项（滑行不变）', main.str(14))


if __name__ == '__main__':
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
