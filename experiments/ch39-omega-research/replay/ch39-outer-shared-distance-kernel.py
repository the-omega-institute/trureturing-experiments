"""两枚留出壳上任意非负、非增的共同距离核的引力弱读数禁区。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


CENTER = 2_306_888_100_085


def weak_bands(low_weak, high_weak, p0, p1):
    # 高截止包围的误差随 x 单调减小。低端 1≤x/p0≤256 用角矩级数覆盖。
    large = high_weak(256*p0, p0, p1)
    whole = large
    for e in range(8):
        step = arb(2)**(e-6)
        for j in range(64):
            ratio = ((64+j)*step).union((65+j)*step)
            whole = whole.union(low_weak(ratio))
    assert arb('-.07479') < whole < arb('-.04764')
    return whole, large


def prefix_cap(row):
    total = arb(0)
    cap = arb(0)
    for k, (_, y) in enumerate(row, 1):
        total += y
        cap = max(cap, (total/k).upper())
    return cap


def certify(number, rows, weights, field, data):
    floor, p0, p1, low_weak, high_weak, whole, large, constants, moments, heat = data
    assert len(rows) == 19 and all(len(row) == 171 for row in rows)
    assert all(c > 0 or c < 0 for c in weights) and weights[15] > 0 > weights[17]
    if_high = sum(c*(large if j == 15 else whole) for j, c in enumerate(weights))
    assert if_high < 0

    sixteenth, eighteenth = rows[15], rows[17]
    assert eighteenth[2][0] < sixteenth[0][0] < sixteenth[2][0] < eighteenth[3][0]
    third = sixteenth[2][1]
    assert third > arb('9.99e-10')
    ratio = 171*16*floor/third
    front = min(y.lower() for _, y in eighteenth[:3])
    lower = front/(1+168*ratio/3)
    assert lower > arb('.75')
    caps = [prefix_cap(row) for row in rows]
    assert all(cap < 1 for cap in caps)
    k, w, error_k, error_w = constants(field)
    variation = sum(abs(c) for c in weights)
    def moment_upper(degree):
        return sum(c*cap**degree for c, cap in zip(weights, caps, strict=True) if c > 0) + weights[17]*lower**degree
    k_lower = -k*moment_upper(2) - variation*error_k
    w_upper = w*moment_upper(4) + variation*error_w
    assert k_lower > 1 and w_upper < 0
    print('编号、第16层高时弱差上界、低时共同核权比与第18层振幅下界、K下界、W上界', number,
          if_high.upper().str(12), ratio.upper().str(10), lower.lower().str(11),
          k_lower.lower().str(12), w_upper.upper().str(12), flush=True)

    # 同一连续族两端对三个目标分别异号，证实每个目标均可单独达到。
    for exponent in (0, 512):
        columns = []
        for row in rows:
            factors = [(row[0][0]/d)**exponent for d, _ in row]
            y = sum(a*v for a, (_, v) in zip(factors, row, strict=True))/sum(factors)
            assert floor < y < 1
            x = y*y
            weak = low_weak(x/p0) if x < 256*p0 else high_weak(x, p0, p1)
            columns.append((*moments(field, heat, x), weak))
        K, W, D = [sum(c*v for c, v in zip(weights, col, strict=True)) for col in zip(*columns, strict=True)]
        if exponent == 0:
            assert K < 1 and W > arb('1e-122') and D < 0
        else:
            assert K > 1 and W < arb('1e-122') and D > 0
        print('同壳逐项匹配的幂指数端点与 K,W,D', exponent, K.str(9), W.str(9), D.str(9), flush=True)


def verify(root):
    training, _ = run_path(str(root/'ch39-three-shell-two-space-readout.py'))['verify'](root, return_data=True)
    native = run_path(str(root/'ch39-native-rh-ray-pvm.py'))
    gravity = run_path(str(root/'ch39-three-index-free-atoms.py'))
    high_weak = run_path(str(root/'ch39-full-neutral-free-cutoff.py'))['weak_column']
    constants = run_path(str(root/'ch39-common-response-required-drawdown.py'))['field_bounds']
    numbers = (CENTER-2, CENTER+2)
    moments, heat, fields = gravity['broken_fields'](root, numbers)
    ctx.prec = 320
    floor = arb('91.1876')/arb('1.22089e19')
    p0, p1 = floor**2, (arb('1.176397e14')/arb('1.22089e19'))**2
    assert len(training) == 171 and all(floor < v < 1 for _, v in training)
    low_weak = native['weak_series']()
    whole, large = weak_bands(low_weak, high_weak, p0, p1)
    for number, field in zip(numbers, fields, strict=True):
        _, levels, pivots, _ = next(x for x in native['shell_weights'](root, number) if x[0] == 6)
        ctx.prec = 320
        rows = []
        for t in levels:
            row = sorted(((abs(t-s), v) for s, v in training), key=lambda item: item[0].mid())
            assert row[0][0] > 0 and all(a[0] < b[0] for a, b in zip(row, row[1:]))
            rows.append(row)
        certify(number, rows, [z.real for z in pivots[0]], field,
                (floor, p0, p1, low_weak, high_weak, whole, large, constants, moments, heat))


if __name__ == '__main__':
    verify(Path(__file__).resolve().parent)
