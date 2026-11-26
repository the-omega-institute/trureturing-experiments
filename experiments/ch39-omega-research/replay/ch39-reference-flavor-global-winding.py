"""精确参考味坐标与最小割三窗同时成立时，连续全局最低点的最小绕数。"""
from collections import Counter
from functools import cache
from pathlib import Path
from runpy import run_path
import json
import sys
from flint import acb, arb, arb_mat, ctx

def lower(root, lib, tri, delta, D):
    """各行正余量与全局极小、两偏导为零不相容。浮点中点只选择有理竞争证书的凸组合。"""
    previous = run_path(str(root / 'ch39-twentyseven-coefficient-window.py'))
    inertia = run_path(str(root / 'ch39-neutral-mean-phase-inertia.py'))
    model = run_path(str(root / 'ch39-mincut-three-window-angle.py'))
    a, b = (arb('18.94088212'), arb('18.94414547'))
    lo, hi = model['sign']['density']['eigen_band']()
    assert a < lo < hi < b
    data = json.loads(Path(__file__).with_name('ch39-reference-flavor-competitors.json').read_text(), parse_float=str)
    coeffs = [lib['coeff'](tri, q, delta) for q in range(1, 36)]
    counts = [0, 0]
    for sign in (-1, 1):
        probes = []
        for origin, (u, v), points in data:
            u, v = (sign * arb(u), sign * arb(v))
            base = []
            for q, c in enumerate(coeffs, 1):
                cost = sum((arb(w) * lib['energy'](q, c, sign * arb(x), sign * arb(y)) for x, y, w in points))
                assert all((arb(w) > 0 for x, y, w in points))
                base.append(cost)
            probes.append((float(origin), u, v, sum((arb(w) for x, y, w in points)), base))
        pending = [(-arb.pi() / 16, arb.pi() / 16)]
        while pending:
            left, right = pending.pop()
            t = left.union(right)
            x = sign * t
            y = sign * (t - D)
            if abs(y) > arb.pi() / 16 or previous['signs'](inertia, 2, a, x, y)[0] >= 16 or previous['signs'](inertia, 2, b, x, y)[1] < 16:
                counts[0] += 1
                continue
            found = False
            energy = [lib['energy'](q, c, x, y) for q, c in enumerate(coeffs, 1)]
            gradients = [lib['derivatives']([(q, c)], x, y)[:2] for q, c in enumerate(coeffs, 1)]
            nearby = sorted(probes, key=lambda z: abs(float(t) - z[0]))
            left_samples = [v for v in probes if v[0] <= float(t)]
            right_samples = [v for v in probes if v[0] > float(t)]
            if left_samples and right_samples:
                l, r = (left_samples[-1], right_samples[0])
                w = arb((float(t) - l[0]) / (r[0] - l[0]))
                assert 0 <= w <= 1
                mixed = (0, *[(1 - w) * a + w * b for a, b in zip(l[1:4], r[1:4])], [(1 - w) * a + w * b for a, b in zip(l[4], r[4])])
                nearby.insert(0, mixed)
            for center, u, v, mass, base in nearby:
                if all((mass * en - ref + u * gx + v * gy > arb('.0001') for en, ref, (gx, gy) in zip(energy, base, gradients))):
                    found = True
                    break
            if found:
                counts[1] += 1
                continue
            mid = (left + right) / 2
            assert right - left > arb('1e-9'), (x, y, center)
            pending.extend(((left, mid), (mid, right)))
    print('完成35圈必要性排除', counts, flush=True)

def witness(lib, tri, delta, d):
    """固定三枚有理权。在全部指定高度上用同一三阶线性系统解剩余三权，并覆盖全环面。"""
    orders = (1, 2, 30, 32, 33, 36)
    coeffs = {q: lib['coeff'](tri, q, delta) for q in orders}
    x = arb('-.01046')
    y = x - d
    locked = {1: arb('.03744'), 2: arb('.02331'), 30: arb('.08198')}
    gradients = {q: lib['derivatives']([(q, c)], x, y)[:2] for q, c in coeffs.items()}
    rhs = [-sum((w * gradients[q][j] for q, w in locked.items())) for j in (0, 1)] + [1 - sum(locked.values())]
    a = arb_mat([[*gradients[q], 1] for q in orders[3:]]).transpose()
    weights = {**locked, **dict(zip(orders[3:], a.solve(arb_mat([[v] for v in rhs])).entries()))}
    assert all((w > arb('.02') for w in weights.values()))
    modes = [(q, [weights[q] * z for z in coeffs[q]]) for q in orders]

    drift = global_cover(lib, tri, weights, modes, x, y)
    return x, y, drift


def global_cover(lib, tri, weights, modes, x, y):
    """用三组一维谐波在整个环面同时包围能量与梯度。"""
    trial = sum(lib["energy"](q, coeffs, x, y) for q, coeffs in modes)
    box = arb(".00002")

    @cache
    def harmonic(axis, n, index):
        if axis == 2:
            t = arb.pi()*arb(2*index, 2)/n
        else:
            t = arb.pi()*(arb(2*index+1, 1)/n-1)
        value, slope = arb(0), arb(0)
        for q, coefficients in modes:
            c = coefficients[axis]
            rotated = c*acb(0, q*t).exp()
            value += c.real-rotated.real
            slope += q*rotated.imag
        return value, slope

    pending = [(1, 0, 0)]
    counts = [0, 0, 0]
    while pending:
        n, i, j = pending.pop()
        ax, by, diag = harmonic(0, n, i), harmonic(1, n, j), harmonic(2, n, i-j)
        if ax[0]+by[0]+diag[0] > trial:
            counts[0] += 1
            continue
        a, b = [arb.pi()*(arb(2*k+1, 1)/n-1) for k in (i, j)]
        if any(abs(a-sign*x)<box and abs(b-sign*y)<box for sign in (-1, 1)):
            counts[1] += 1
            continue
        if not (ax[1]+diag[1]).contains(0) or not (by[1]-diag[1]).contains(0):
            counts[2] += 1
            continue
        assert n < 2**23, (n, i, j, a, b)
        pending.extend((2*n, 2*i+u, 2*j+v) for u in (0, 1) for v in (0, 1))
    for t in (-1, 1):
        a, b = (t * x + arb(0, box), t * y + arb(0, box))
        _, _, h, j, k = lib['derivatives'](modes, a, b)
        assert h > 50 and (h - 50) * (k - 50) > j * j
    drift = arb(0)
    for i, j in ((1, 0), (0, 1), (1, -1)):
        angle = i * x + j * y + arb(0, box * (abs(i) + abs(j)))
        for q, w in weights.items():
            drift -= 2 * w * sum(((z ** q).imag for z, p in tri[i, j] if p)) * (q * angle).sin()
    print('正权项数', len(weights), '能量、极小点方框、非驻点排除格', counts,
          '所选支减去镜像的能量差除以辐边亏额', drift.str(12), flush=True)
    return drift


def conjugation_symmetry(root):
    """单位辐边的谐波均有精确的镜像对。"""
    faces = run_path(str(root / 'ch39-mincut-real-triangle-tomography.py'))
    table = faces['enumerate_faces']()
    for frequency in ((1, 0), (0, 1), (1, -1)):
        counts = Counter((z for k, z in table.values() if k == frequency))
        assert counts == Counter({faces['conjugate'](z): n for z, n in counts.items()})


def physical_all_heights(root, x, y, cap):
    physical = run_path(str(root / "ch39-first-three-mode-observables.py"))
    model, initial = physical["model_inputs"]()
    chiral = physical["joint"]["chiral"]
    phases = {724: acb(0, x).exp(), 804: acb(0, y).exp()}
    rows = chiral["source"]["spectral_supplier"](model)(arb("1e13"), clocks=phases)
    # 任意两种径向刚度的三辐边差为带相位的三星图，算子范数至多 4δ_max。
    # 拉普拉斯的谱至多 36（19 站、最大度至多 18），中性中心至多 15。
    # v|v-c| 的 Lipschitz 界为 (2*36+15)|Δv|+36|Δc|。
    err = (2*36+15)*4*cap + 36*cap/3
    assert err < arb("2e-6")
    bands = {q: chiral["ordered"]([v+arb(0, err) for v, charge in rows if charge == q])
             for q in (0, *model["charges"])}
    assert all(bands[q][2] < 9 and bands[q][3] > 10 for q in model["charges"])
    _, alpha, curve = physical["layer_observables"](model, initial, bands)
    assert arb(".117") < alpha < arb(".119") and arb(".998") < curve < arb("1.002")
    print("全高度的强耦合与曲率比", alpha.str(14), curve.str(14))


def verify(root):
    lib = run_path(str(root / 'ch39-mincut-global-winding.py'))
    shell = run_path(str(root / 'ch39-three-ray-chiral-boundary.py'))
    tri = lib['terms'](shell)
    ctx.prec = 224
    conjugation_symmetry(root)
    cap = arb(23)/25*arb(2)**-20/14**2
    delta = arb(0).union(cap)
    phases = shell['geo']['phases'](acb.zeta_zero(11661).imag)
    d = (phases[724] / phases[804]).arg()
    flavor = run_path(str(root / 'ch39-chiral-neutral-flavor-phase.py'))
    a, b = [flavor['singlet_at']({724: z, 804: acb(1)})[1] for z in (acb(1), acb(0, 1))]
    alternative = 2 * (b / a).atan() - d
    assert a > 0 and abs(b) < arb('1e-45') and (acb.zeta_zero(1).imag > 14)
    both = abs(d).union(abs(alternative))
    full = arb(0).union(arb('4.48e-9'))
    assert arb(23) / 25 * arb(2) ** (-20) / 14 ** 2 < arb('4.48e-9') and shell['soft']['point']['study']['offset'] > 0
    lower(root, lib, tri, full, both)
    x, y, drift = witness(lib, tri, delta, alternative)
    assert drift < -3.6
    physical_all_heights(root, x, y, cap)
    ref = flavor['singlet'](11661)[1]
    found = flavor['singlet_at']({724: acb(0, x).exp(), 804: acb(0, y).exp()})[1]
    assert found.overlaps(ref)
    print('参考有向角距、精确等味线上的相位与坐标', d.str(15), x.str(15), y.str(15), found.str(15))
if __name__ == '__main__':
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    verify(root)
