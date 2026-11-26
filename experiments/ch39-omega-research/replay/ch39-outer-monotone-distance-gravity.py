"""三壳171项校准数据的非增距离加权，在两处留出壳分别构造三项引力弱匹配。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


CENTER = 2_306_888_100_085


def inner_prefix_range(t, training):
    ordered = sorted(((abs(t - node), y) for node, y in training), key=lambda row: row[0].mid())
    assert ordered[0][0] > 0
    assert all(left[0] < right[0] for left, right in zip(ordered, ordered[1:]))
    total = arb(0)
    means = []
    for k, (_, y) in enumerate(ordered, 1):
        total += y
        means.append(total / k)
    # 此处需要可达内区间。下端取某个均值的上界，上端取某个均值的下界。
    return min(z.upper() for z in means), max(z.lower() for z in means)


def verify_outer(number, weights, ranges, field, moments, heat, low_weak, high_weak):
    mz, mp = arb('91.1876'), arb('1.22089e19')
    p0, p1 = (mz/mp)**2, (arb('1.176397e14')/mp)**2
    middle = arb('.15')
    low_weight = weights[0] + weights[15]
    mid_weight = sum(weights[1:15]) + weights[16]
    tail_weights = weights[17:19]
    assert len(weights) == 19 and sum(weights).contains(1)
    assert low_weight > 0 and tail_weights[0] < 0 < tail_weights[1]
    t_lo, t_hi = arb('1.25'), arb(3)
    for j, (lo, hi) in enumerate(ranges[:17]):
        if j in (0, 15):
            assert lo < (p0*t_lo).sqrt() and (p0*t_hi).sqrt() < hi
        else:
            assert lo < middle < hi

    mid_kw = moments(field, heat, middle**2)
    target = [arb(1), arb('1e-122')]
    def residual(point):
        tails = [moments(field, heat, x) for x in point]
        return arb_mat([[mid_weight*mid_kw[i] + sum(w*v[i] for w, v in zip(tail_weights, tails, strict=True)) - goal]
                        for i, goal in enumerate(target)])

    def jacobian(point):
        cols = []
        for weight, x in zip(tail_weights, point, strict=True):
            decay = [(-z/x).exp() for _, _, z in field]
            k = sum(a*e for (_, a, _), e in zip(field, decay, strict=True))/(2*arb.pi())
            w = -x*sum(g*e for (g, _, _), e in zip(field, decay, strict=True))/(4*arb.pi())
            cols.append([weight*k, weight*w])
        return arb_mat(cols).transpose()

    point = [arb('.12'), arb('.45')]
    for _ in range(6):
        step = jacobian(point).solve(residual(point)).entries()
        point = [(x-dx).mid() for x, dx in zip(point, step, strict=True)]
    radius = arb(2)**-35
    box = [x + arb(0, radius) for x in point]
    for (lo, hi), x in zip(ranges[17:], box, strict=True):
        assert lo < x.sqrt() < hi
    # 对整个低截止区间应用 0<E₂≤1 与 0<E₃≤1/2 的积分界。
    assert all(z >= 0 for _, _, z in field)
    max_low = p0*t_hi
    low_error = [abs(low_weight)*max_low*sum(abs(a) for _, a, _ in field)/(2*arb.pi()),
                 abs(low_weight)*max_low**2*sum(abs(g) for g, _, _ in field)/(8*arb.pi())]
    remainder = arb_mat([[arb(0, v.upper())] for v in low_error])
    pre = jacobian(point).inv()
    assert not pre.det().contains(0)
    defect = arb_mat([[1, 0], [0, 1]]) - pre*jacobian(box)
    shift = pre*(residual(point) + remainder)
    ratio = max((abs(s)/radius + sum(abs(z) for z in row)).upper()
                for s, row in zip(shift.entries(), defect.tolist(), strict=True))
    assert ratio < 1

    fixed = mid_weight*high_weak(middle**2, p0, p1)
    fixed += sum(w*high_weak(x, p0, p1) for w, x in zip(tail_weights, box, strict=True))
    def weak_face(t):
        return low_weight*low_weak(t) + fixed
    assert weak_face(t_lo) > 0 > weak_face(t_hi)
    for _ in range(25):
        mid = (t_lo+t_hi)/2
        value = weak_face(mid)
        if value > 0:
            t_lo = mid
        elif value < 0:
            t_hi = mid
        else:
            break
    d_lo, d_hi = weak_face(t_lo), weak_face(t_hi)
    assert d_lo > 0 > d_hi
    print('编号、共同低截止/GeV、末两层平方截止、热核面相对误差、弱差两面', number,
          (mz*t_lo.union(t_hi).sqrt()).str(12), *[v.str(12) for v in box], ratio.str(7),
          d_lo.str(7), d_hi.str(7), flush=True)
    low = (p0*t_lo.union(t_hi)).sqrt()
    return [low if j in (0, 15) else middle for j in range(17)] + [x.sqrt() for x in box]


def verify(root, *, training=None):
    if training is None:
        training, _ = run_path(str(root/'ch39-three-shell-two-space-readout.py'))['verify'](root, return_data=True)
    assert len(training) == 171
    native = run_path(str(root/'ch39-native-rh-ray-pvm.py'))
    gravity = run_path(str(root/'ch39-three-index-free-atoms.py'))
    high_weak = run_path(str(root/'ch39-full-neutral-free-cutoff.py'))['weak_column']
    numbers = (CENTER-2, CENTER+2)
    moments, heat, fields = gravity['broken_fields'](root, numbers)
    ctx.prec = 320
    weak = native['weak_series']()
    responses = {}
    for number, field in zip(numbers, fields, strict=True):
        _, levels, pivots, _ = next(x for x in native['shell_weights'](root, number) if x[0] == 6)
        ctx.prec = 320
        ranges = [inner_prefix_range(t, training) for t in levels]
        responses[number] = verify_outer(number, [w.real for w in pivots[0]], ranges, field, moments, heat, weak, high_weak)
    return responses


if __name__ == '__main__':
    verify(Path(__file__).resolve().parent)
