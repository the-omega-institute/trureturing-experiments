"""实际相邻壳任意站位对角联络的弱差阻碍，以及两角同时重建两项引力矩。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


CENTER = 2_306_888_100_085


def pole_bound(source, target, src, dest):
    terms = [abs(target[i, dest])*abs(source[i, src]) for i in range(19)]
    return terms[14] - sum(v for i, v in enumerate(terms) if i != 14)


def two_angle_solution(source, target, amplitudes, weights, field, moments, heat, seed, radius):
    full = target.conjugate().transpose()*source
    rows = []
    for j in range(19):
        row = []
        for i in range(19):
            a, b = [target[s,j].conjugate()*source[s,i] for s in (6,14)]
            row.append((full[j,i]-a-b, a, b))
        rows.append(row)
    base = min(v.lower() for v in amplitudes)
    delta = [a-base for a in amplitudes]
    def value_and_jacobian(angles, *, values=True):
        phases = [acb(0, x).exp() for x in angles]
        residual = [-arb(1), -arb('1e-122')]
        jac = arb_mat(2, 2)
        for row, weight in zip(rows, weights, strict=True):
            f, df, dg = base, arb(0), arb(0)
            for (z, a, b), dy in zip(row, delta, strict=True):
                p, q = phases[0]*a, phases[1]*b
                total = z+p+q
                f += dy*abs(total)**2
                df -= 2*dy*(total.conjugate()*p).imag
                dg -= 2*dy*(total.conjugate()*q).imag
            assert f > 0
            x = f**2
            if values:
                k, w = moments(field, heat, x)
                residual[0] += weight*k
                residual[1] += weight*w
            decay = [(-z/x).exp() for _, _, z in field]
            dk = sum(a*e for (_, a, _), e in zip(field, decay, strict=True))/(2*arb.pi())
            dw = -x*sum(g*e for (g, _, _), e in zip(field, decay, strict=True))/(4*arb.pi())
            for j, deriv in enumerate((df, dg)):
                jac[0,j] += 2*weight*f*deriv*dk
                jac[1,j] += 2*weight*f*deriv*dw
        return arb_mat([[v] for v in residual]), jac
    point = list(map(arb, seed))
    for _ in range(4):
        res, jac = value_and_jacobian(point)
        point = [(x-dx).mid() for x, dx in zip(point, jac.solve(res).entries(), strict=True)]
    radius = arb(radius)
    box = [v+arb(0,radius.upper()) for v in point]
    res, jac = value_and_jacobian(point)
    pre = jac.inv()
    assert not pre.det().contains(0)
    _, wide = value_and_jacobian(box, values=False)
    error = arb_mat([[1,0],[0,1]])-pre*wide
    ratio = max((abs(z)/radius+sum(abs(e) for e in row)).upper()
                for z, row in zip((pre*res).entries(), error.tolist(), strict=True))
    assert ratio < 1
    print('两项热核的站6、14相位矩形与自映射相对位移',
          *[x.mid().str(14) for x in point], radius.str(4), ratio.str(8), flush=True)
    phases = [acb(0, x).exp() for x in box]
    return [base + sum(dy*abs(z+phases[0]*a+phases[1]*b)**2
                       for (z, a, b), dy in zip(row, delta, strict=True)) for row in rows]


def verify(root):
    # 两侧荷六训练值由原证书直接给出，首层左端采用两次转向的构造。
    left = run_path(str(root/'ch39-left-cutoff-two-turns.py'))['verify'](root)
    right = run_path(str(root/'ch39-adjacent-ordered-charge-six-cutoffs.py'))['verify'](root)[2]
    shared = run_path(str(root/'ch39-adjacent-pvm-transport.py'))
    shell = run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    pivot = run_path(str(root/'ch39-native-rh-ray-pvm.py'))
    bounds = run_path(str(root/'ch39-outer-shared-distance-kernel.py'))['weak_bands']
    high_weak = run_path(str(root/'ch39-full-neutral-free-cutoff.py'))['weak_column']
    gravity = run_path(str(root/'ch39-three-index-free-atoms.py'))
    numbers = [CENTER+i for i in (-2, 2, -3, 3)]
    moments, heat, fields = gravity['broken_fields'](root, numbers)
    by_number = dict(zip(numbers, fields, strict=True))
    ctx.prec = 320
    mz_ratio = arb('91.1876')/arb('1.22089e19')
    whole, large = bounds(pivot['weak_series'](), high_weak, mz_ratio**2, (arb('1.176397e14')/arb('1.22089e19'))**2)
    for side, source_data, seed in zip((-1,1), (left,right),
            (('-.680146','.506372'),('-.692332','-.753392')), strict=True):
        origin, dest = CENTER+side, CENTER+2*side
        source_levels, source = shared['spectrum'](shell, 6, origin)
        _, target = shared['spectrum'](shell, 6, dest)
        assert all(a.overlaps(b) for a, b in zip(source_levels, source_data[0], strict=True))
        amplitudes = [x.sqrt() for x in source_data[1]]
        assert all(mz_ratio < a < 1 for a in amplitudes)
        _, _, cross, _ = next(t for t in pivot['shell_weights'](root, dest) if t[0] == 6)
        ctx.prec = 320
        weights = [w.real for w in cross[0]]
        c = pole_bound(source, target, 17, 15)
        assert c > 0
        min_f = amplitudes[17]*c**2
        distance = c - (16*mz_ratio/amplitudes[17]).sqrt()
        assert min_f > 16*mz_ratio and distance > (arb('.20364') if side < 0 else arb('.09590'))
        weak_bound = sum(w*(large if j == 15 else whole) for j,w in enumerate(weights))
        assert weak_bound < 0
        print('编号、18到16最小重叠模、第16层振幅下界、弱差上界、非对角算子距离门槛', dest,
              c.str(13), min_f.str(12), weak_bound.upper().str(12), distance.str(13), flush=True)
        first = two_angle_solution(source, target, amplitudes, weights, by_number[dest], moments, heat, seed, '2e-9')
        final = CENTER+3*side
        _, endpoint = shared['spectrum'](shell, 6, final)
        _, _, last_cross, _ = next(t for t in pivot['shell_weights'](root, final) if t[0] == 6)
        ctx.prec = 320
        last_weights = [w.real for w in last_cross[0]]
        # 各步均取站位14的反三角界。两条中继路径可以直接相加，因为每一步输运谱振幅。
        links = [(max(arb(0), pole_bound(source, target, 17, j).lower()),
                  max(arb(0), pole_bound(target, endpoint, j, 15).lower())) for j in (15, 17)]
        min_end = amplitudes[17]*sum(h1**2*h2**2 for h1, h2 in links)
        end_weak = sum(w*(large if j == 15 else whole) for j, w in enumerate(last_weights))
        assert min_end > 16*mz_ratio and end_weak < 0
        print('两步编号、中继第16、18层的两段模下界、末端第16层振幅下界与弱差上界', final,
              [[z.str(9) for z in path] for path in links], min_end.lower().str(11), end_weak.upper().str(11), flush=True)
        seed2 = ('-.034791','.418339') if side < 0 else ('-.272280','.811764')
        two_angle_solution(target, endpoint, first, last_weights, by_number[final], moments, heat, seed2, '1e-5')


if __name__ == '__main__':
    verify(Path(__file__).resolve().parent)
