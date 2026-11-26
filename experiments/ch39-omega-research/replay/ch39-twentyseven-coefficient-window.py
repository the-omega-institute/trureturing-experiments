"""第一与第二十七圈的非负比例、单位辐边及高处软壳的全局真空窗。"""
from collections import Counter
from pathlib import Path
from runpy import run_path
import sys
from flint import acb, acb_mat, arb, arb_mat, ctx

def outer(root, lib, coeffs, lo, hi):
    """先用荷二与荷六惯性限制物理区域，再联立归一混合的两个驻点条件和全局能量比较。"""
    model = run_path(str(root / 'ch39-mincut-three-window-angle.py'))
    ctx.prec = 176
    band = model['sign']['density']['eigen_band']()
    lower, upper = (arb('18.94088212'), arb('18.94414547'))
    assert lower < band[0] < band[1] < upper
    c, d = [model['cell_count'](2, v) for v in (lower, upper)]
    six = model['cell_count'](6, (15 - arb(185).sqrt()) / 2 + arb('3e-8'))
    cnt = 0
    pending = [(1, 0, 0)]
    while pending:
        n, i, j = pending.pop()
        if c(n, i, j)[0] >= 16 or d(n, i, j)[1] < 16 or six(n, i, j)[1] < 1:
            continue
        x, y = [arb.pi() / 16 * (arb(2 * k + 1, 1) / n - 1) for k in (i, j)]
        g, f = [lib['derivatives']([(q, z)], x, y)[:2] for q, z in zip((1, 27), coeffs)]
        a = arb(0).union(arb(1))
        for b, v in zip(g, f):
            if b > 0 and v > 0 or (b < 0 and v < 0):
                a = None
                break
            if not (b - v).contains(0):
                cut = b / (b - v)
                if not a.overlaps(cut):
                    a = None
                    break
                a = a.intersection(cut)
        if a is not None:
            xx, yy = (arb('.071'), arb('.076'))
            d1, d2 = [lib['energy'](q, z, x, y) - lib['energy'](q, z, xx, yy) for q, z in zip((1, 27), coeffs)]
            if (1 - a) * d1 + a * d2 > 0:
                a = None
        if a is None or lo / (1 + lo) < a < hi / (1 + hi):
            cnt += 1
            continue
        assert n < 2 ** 22, (n, i, j, a, a / (1 - a))
        pending.extend(((n * 2, i * 2 + u, j * 2 + v) for u in (0, 1) for v in (0, 1)))
    print('全局合格比例的必要开区间、驻点或能量终止格', lo, hi, cnt)

def cover(lib, coeffs, r):
    """用整个比例段上的能量与梯度筛去全环面格子，余区的两对角项和钟差项给出严格正定 Hessian。"""
    xx, yy = (arb('.071'), arb('.076'))
    ref = [lib['energy'](q, c, xx, yy) for q, c in zip((1, 27), coeffs)]
    modes = [(1, coeffs[0]), (27, [v * r for v in coeffs[1]])]
    pending = [(1, 0, 0)]
    while pending:
        n, i, j = pending.pop()
        x, y = [arb.pi() * (arb(2 * k + 1, 1) / n - 1) for k in (i, j)]
        if any((arb('.065') < s * x < arb('.077') and arb('.071') < s * y < arb('.082') for s in (-1, 1))):
            continue
        dx, dy = [lib['energy'](q, c, x, y) - t for q, c, t in zip((1, 27), coeffs, ref)]
        if dx + r * dy > 0:
            continue
        gx, gy, *_ = lib['derivatives'](modes, x, y)
        if not gx.contains(0) or not gy.contains(0):
            continue
        assert n < 2 ** 18, (n, i, j, x, y)
        pending.extend(((n * 2, i * 2 + u, j * 2 + v) for u in (0, 1) for v in (0, 1)))
    for s in (-1, 1):
        x = s * arb('.065').union(arb('.077'))
        y = s * arb('.071').union(arb('.082'))
        diagonal = [sum((q * q * (c[k] * acb(0, q * t).exp()).real for q, c in modes)) for k, t in enumerate((x, y, x - y))]
        assert diagonal[0] > 5 and diagonal[1] > 5 and (diagonal[2] > 0)
    print('全环面只余两个严格凸方框', flush=True)

def enclose(lib, coeffs, r, sign):
    """对各个固定比例，预条件梯度映射将相位方框严格映回自身，凸性保证其零点惟一。"""

    def deriv(p, t, x, y):
        g, h, a, b, d = lib['derivatives']([(p, t)], x, y)
        return (arb_mat([[g], [h]]), arb_mat([[a, b], [b, d]]))

    def fixed(a, b, R):
        g, H = deriv(1, coeffs[0], a, b)
        f, J = deriv(27, coeffs[1], a, b)
        return (g + R * f, H + R * J)
    x, y = (sign * arb('.071'), sign * arb('.076'))
    for _ in range(6):
        g, H = fixed(x, y, r.mid())
        step = H.solve(g)
        x = (x - step[0, 0]).mid()
        y = (y - step[1, 0]).mid()
    g, H = deriv(1, coeffs[0], x, y)
    f, J = deriv(27, coeffs[1], x, y)
    C = (H + r.mid() * J).inv()
    width = (r.upper() - r.lower()) * arb('.22') + arb('1e-13')
    xx, yy = (x + arb(0, width), y + arb(0, width))
    _, actual = fixed(xx, yy, r)
    residual = C * g + r * (C * f)
    remainder = arb_mat([[1, 0], [0, 1]]) - C * actual
    for i in (0, 1):
        assert abs(residual[i, 0]) + width * sum((abs(remainder[i, j]) for j in (0, 1))) < width
    return (xx, yy)

def reference(shell, q, x, y):
    """计算用于选择谱分隔点的参考谱。分隔点的接受依赖整格惯性证书。"""
    mat = acb_mat(19, 19)
    phases = {724: acb(0, q * x).exp(), 804: acb(0, q * y).exp()}
    for i, j, d, z in shell['geo']['edges']:
        v = z ** q * phases.get(d, 1)
        mat[i, i] += 1
        mat[j, j] += 1
        mat[i, j] = -v
        mat[j, i] = -v.conjugate()
    return sorted((z.real for z in mat.eig(multiple=True)), key=float)

def signs(inertia, q, level, x, y):
    """三条轮换特征标各用中心一阶展开和二阶余项包围，给出整格上阈下本征值个数的上下界。"""
    count, cs = inertia['schur_coefficients'](q, level)
    centers = (x.mid(), y.mid(), y.mid() - x.mid())
    phases = [acb(0, q * t).exp() for t in centers]
    vals = []
    for c, a, b, d in cs:
        modes = [z * p for z, p in zip((a, b, d), phases)]
        za, zb, zd = modes
        mid = c + 2 * (za + zb + zd).real
        rx, ry = (q * x.rad(), q * y.rad())
        err = 2 * (abs((za - zd).imag) * rx + abs((zb + zd).imag) * ry) + abs(a) * rx ** 2 + abs(b) * ry ** 2 + abs(d) * (rx + ry) ** 2
        vals.append(mid + arb(0, err.abs_upper()))
    return (16 - count + sum((v < 0 for v in vals)), 16 - count + sum((not v > 0 for v in vals)))

def bands(shell, inertia, q, x, y, delta):
    """对每个有序本征值认证上下两个分隔点，再按单位辐边的谱均值 15 折叠并重新排序。"""
    z = reference(shell, q, x.mid(), y.mid())
    corners = [reference(shell, q, p1, p2) for p1 in (x.lower(), x.upper()) for p2 in (y.lower(), y.upper())]
    vals = []
    for j, v in enumerate(z):
        ends = []
        for sg in (-1, 1):
            gap = max([arb('1e-10'), *[sg * (ww[j] - v) for ww in corners]])
            for _ in range(10):
                level = v.mid() + sg * gap * arb('1.5')
                lo, hi = signs(inertia, q, level, x, y)
                if sg < 0 and hi <= j or (sg > 0 and lo > j):
                    ends.append(level)
                    break
                gap *= 2
            else:
                raise RuntimeError((q, j, sg, gap))
        vals.append(ends[0].union(ends[1]))
    # 三辐边星图的最大本征值为 4；软化使每个有序能级下降至多 4δ。
    center = 15-delta/3
    shifted = [v-4*delta for v in vals]
    masses = [v*abs(v-center) for v in shifted]
    return [a.union(b) for a, b in zip(sorted((z.lower() for z in masses)), sorted((z.upper() for z in masses)), strict=True)]

def flavor(root, lib, coeff, tri):
    """在整个必要比例段上控制两个全局极小分支的相对味读数。"""
    cp = run_path(str(root / 'ch39-three-window-two-probability-cp.py'))['cp']
    ctx.prec, ctx.cap = (384, 3)
    f = cp['flavor']
    data = cp['two']['spectral_data']()
    polys = cp['cp_data'](data, f['reference_matrix']())
    ref = f['singlet'](11661)[1]
    cu = cp['two']['probabilities'](ref, data)[0][0]
    cj = cp['cp_value'](ref, polys)
    ax, bx = [f['singlet_at']({724: z, 804: acb(1)})[1] for z in (acb(1), acb(0, 1))]
    # 精确几何系数的共轭配对给出 δ=0 时的镜面对称。
    faces = run_path(str(root/'ch39-mincut-real-triangle-tomography.py'))
    face_table = faces['enumerate_faces']()
    soft = {}
    for k in ((1,0),(0,1),(1,-1)):
        counts = Counter(z for mode,z in face_table.values() if mode==k)
        assert counts == Counter({faces['conjugate'](z):n for z,n in counts.items()})
        soft[k] = [sum((z**q).imag for z,n in tri[k] if n==1) for q in (1,27)]
    for sign in (-1, 1):
        upp = []
        jrr = []
        for i in range(1240, 1410):
            x, y = enclose(lib, coeff, (arb(i) / 10000).union(arb(i + 1) / 10000), sign)
            assert arb('.065') < sign * x < arb('.077') and arb('.071') < sign * y < arb('.082')
            r = (arb(i)/10000).union(arb(i+1)/10000)
            split = arb(0)
            for k,t in zip(((1,0),(0,1),(1,-1)),(sign*x,sign*y,sign*(x-y)),strict=True):
                one,many = soft[k]
                split -= 2*(one*t.sin()+r*many*(27*t).sin())
            assert split < arb("-.04")
            d = x - y
            c = ax * d.cos() + bx * d.sin()
            uc = cp['two']['probabilities'](c, data)[0]
            uj = cp['cp_value'](c, polys) / cj
            upp.append(uc[0] / cu)
            jrr.append(uj)
        assert all((arb('.2690917') < u < arb('.2692115') for u in upp))
        assert all((arb('.3133773') < j < arb('.3136154') for j in jrr))
        print('相位分支及 U、J 相对参考值的区间', sign, *[(min((t.lower() for t in arr)).str(15), max((t.upper() for t in arr)).str(15)) for arr in [upp, jrr]])

def verify(root, lib):
    shell = run_path(str(root / 'ch39-three-ray-chiral-boundary.py'))
    physical = run_path(str(root / 'ch39-first-three-mode-observables.py'))
    inertia = run_path(str(root / 'ch39-neutral-mean-phase-inertia.py'))
    inputs, initial = physical['model_inputs']()
    ctx.prec = 224
    tri = lib['terms'](shell)
    assert shell['soft']['point']['study']['offset'] > 0
    delta = arb(0).union(arb(23)/25*arb(2)**-20/10**24)
    coeffs = [lib['coeff'](tri,q,delta) for q in (1,27)]
    outer(root, lib, coeffs, arb('.124'), arb('.141'))
    interval = arb('.124').union(arb('.141'))
    cover(lib, coeffs, interval)
    flavor(root, lib, coeffs, tri)
    # 已有中性高度谱的 1e-12 余量覆盖单位辐边极限，且不依赖两钟相位。
    chiral = physical['joint']['chiral']
    supply = chiral['source']['spectral_supplier'](inputs)
    neutral = chiral['ordered']([v for v, q in supply(arb('1e13'), clocks={724: acb(1), 804: acb(1)}) if q == 0])
    alphas = []
    curves = []
    for i in range(125, 140):
        for s in (-1, 1):
            for lo, hi in [(2 * i, 2 * i + 1), (2 * i + 1, 2 * i + 2)] if i in (125, 139) else [(2 * i, 2 * i + 2)]:
                r = (arb(lo) / 2000).union(arb(hi) / 2000)
                x, y = enclose(lib, coeffs, r, s)
                assert arb('.065') < s * x < arb('.077') and arb('.071') < s * y < arb('.082')
                data = {q: bands(shell, inertia, q, x, y, delta) for q in inputs['charges']}
                data[0] = neutral
                assert all((data[q][2] < 9 and data[q][3] > 10 for q in inputs['charges']))
                _, alpha, curve = physical['layer_observables'](inputs, initial, data)
                assert arb('.117') < alpha < arb('.119') and arb('.998') < curve < arb('1.002'), (lo, hi, s, alpha, curve)
                alphas += [alpha]
                curves += [curve]
    print('0.125 至 0.140 全闭段上的强耦合与曲率比包围', *[(min((v.lower() for v in seq)).str(13), max((v.upper() for v in seq)).str(13)) for seq in (alphas, curves)])
if __name__ == '__main__':
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    source = Path(sys.argv[2]) if len(sys.argv) > 2 else root / 'ch39-mincut-global-winding.py'
    verify(root, run_path(str(source)))
