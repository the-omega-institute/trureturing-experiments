"""用 Schur–Horn 主序在保全荷六训练振幅总和的两处相邻留出壳实现三项匹配。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


CENTER = 2_306_888_100_085


def ky_fan_gaps(source, output):
    """十九项总和由构造恒等。正向及补集求和任选较紧的严格下界。"""
    n = len(source)
    assert n == len(output) == 19
    source_top = sorted((v.lower() for v in source), reverse=True)
    output_top = sorted((v.upper() for v in output), reverse=True)
    source_bottom = sorted(v.upper() for v in source)
    output_bottom = sorted(v.lower() for v in output)
    gaps = [max((sum(source_top[:k]) - sum(output_top[:k])).lower(),
                (sum(output_bottom[:n-k]) - sum(source_bottom[:n-k])).lower())
            for k in range(1, n)]
    assert min(gaps) > 0, gaps
    return gaps


def certify(number, source, weights, field, moments, heat, low_weak, high_weak):
    pi = arb.pi()
    mz, mp = arb('91.1876'), arb('1.22089e19')
    p0, p1 = (mz / mp)**2, (arb('1.176397e14') / mp)**2
    total = sum(source)
    # 第 1、16 层取同一低振幅，第 18、19 层各留一个热核变量。
    low_weight = weights[0] + weights[15]
    middle_weight = sum(weights[1:15]) + weights[16]
    tail_weights = weights[17:]
    assert (low_weight + middle_weight + sum(tail_weights)).contains(1)
    assert all(z >= 0 for _, _, z in field)

    # 低截止参数约束宽于最终的弱差根。其贡献按积分界全程包围。
    tlo, thi = arb('1.25'), arb(3)
    def low_amp(t):
        return (p0*t).sqrt()
    low = low_amp(tlo.union(thi))
    low_kw = [abs(low_weight)*p0*thi*sum(abs(a) for _, a, _ in field)/(2*pi),
              abs(low_weight)*(p0*thi)**2*sum(abs(g) for g, _, _ in field)/(8*pi)]

    def fill(pair):
        return (total - 2*low - sum(pair)) / 15

    def values(pair):
        midcol = moments(field, heat, fill(pair)**2)
        tails = [moments(field, heat, a**2) for a in pair]
        return arb_mat([[middle_weight*midcol[i] + arb(0, low_kw[i].upper())
                         + sum(w*col[i] for w, col in zip(tail_weights, tails, strict=True)) - goal]
                        for i, goal in enumerate((arb(1), arb('1e-122')))])

    def slope(a):
        decay = [(-z/a**2).exp() for _, _, z in field]
        return (a*sum(v*e for (_, v, _), e in zip(field, decay, strict=True))/pi,
                -a**3*sum(g*e for (g, _, _), e in zip(field, decay, strict=True))/(2*pi))

    def jacobian(pair):
        base = slope(fill(pair))
        return arb_mat([[w*d - middle_weight*b/15 for d, b in zip(slope(a), base, strict=True)]
                        for w, a in zip(tail_weights, pair, strict=True)]).transpose()

    center = [arb('.35'), arb('.67')]
    for _ in range(7):
        delta = jacobian(center).solve(values(center)).entries()
        center = [(a-v).mid() for a, v in zip(center, delta, strict=True)]
    radius = arb('1e-8')
    box = [a + arb(0, radius.upper()) for a in center]
    assert low < fill(box) < box[0] < box[1] < 1
    pre = jacobian(center).inv()
    assert not pre.det().contains(0)
    defect = arb_mat([[1, 0], [0, 1]]) - pre*jacobian(box)
    shift = pre*values(center)
    kratio = max((abs(v)/radius + sum(abs(e) for e in row)).upper()
                 for v, row in zip(shift.entries(), defect.tolist(), strict=True))
    assert kratio < 1, kratio
    fixed = middle_weight*high_weak(fill(box)**2, p0, p1)
    fixed += sum(w*high_weak(a**2, p0, p1) for w, a in zip(tail_weights, box, strict=True))
    def weak_face(t):
        return low_weight*low_weak(t) + fixed
    assert weak_face(tlo) > 0 > weak_face(thi)
    for _ in range(28):
        mid = (tlo + thi)/2
        face = weak_face(mid)
        if face > 0:
            tlo = mid
        elif face < 0:
            thi = mid
        else:
            break
    assert weak_face(tlo) > 0 > weak_face(thi)
    low = low_amp(tlo.union(thi))
    middle = fill(box)
    output = [low if j in (0, 15) else middle for j in range(17)] + box
    assert all(mz/mp < z < 1 for z in output)
    gaps = ky_fan_gaps(source, output)
    print('编号、源振幅总和、最低截止/GeV、中间振幅、末两层振幅、Krawczyk 比、十八个主序差的最小下界',
          number, *[v.str(12) for v in (total, low*mp, middle, *box, kratio, min(gaps))], flush=True)
    print('十八个降序部分和差的严格下界', [v.str(8) for v in gaps], flush=True)
    return output


def verify(root):
    training, _ = run_path(str(root/'ch39-three-shell-two-space-readout.py'))['verify'](root, return_data=True)
    assert len(training) == 171
    gravity = run_path(str(root/'ch39-three-index-free-atoms.py'))
    pivots = run_path(str(root/'ch39-native-rh-ray-pvm.py'))
    transport = run_path(str(root/'ch39-adjacent-pvm-transport.py'))
    weak_high = run_path(str(root/'ch39-full-neutral-free-cutoff.py'))['weak_column']
    targets = (CENTER-2, CENTER+2)
    moments, heat, fields = gravity['broken_fields'](root, targets)
    ctx.prec = 320
    weak_low = pivots['weak_series']()
    outputs = {}
    for target, field, origin in zip(targets, fields, (CENTER-1, CENTER+1), strict=True):
        _, levels, _, _ = next(x for x in pivots['shell_weights'](root, origin) if x[0] == 6)
        source = transport['training_values'](levels, training)
        _, _, cross, _ = next(x for x in pivots['shell_weights'](root, target) if x[0] == 6)
        ctx.prec = 320
        outputs[target] = certify(target, source, [w.real for w in cross[0]], field, moments, heat,
                                  weak_low, weak_high)
    return outputs


if __name__ == '__main__':
    verify(Path(__file__).resolve().parent)
