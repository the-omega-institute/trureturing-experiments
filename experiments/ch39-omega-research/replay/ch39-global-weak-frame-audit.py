"""以精确有理数复核六支舍入与整数角检验，并用投影内积重算两壳弦长。"""
from fractions import Fraction
from importlib import import_module
from unittest.mock import patch

import numpy as np
from flint import acb_mat, arb


model = import_module('ch39-global-independent-weak-frame-separation')


def replay_with_exact_endpoints():
    original = model.coeff
    counts = [0, 0]

    def audited(a, b, c):
        x, y, z = original(a, b, c)
        assert x <= a.lower().fmpq() and y <= b.lower().fmpq()
        assert z >= c.upper().fmpq()
        counts[0] += 2
        counts[1] += 1
        return x, y, z

    with patch.object(model, 'coeff', audited):
        model.verify()
    assert counts[0] == 2 * counts[1] and counts[1] > 0
    print('逐项与精确有理数比较的下端点、上端点', *counts, flush=True)


def chord_distances():
    def projectors(number, levels):
        up, _ = model.raw['native'](model.ROOT, number)
        result = {}
        for level in levels:
            row = up[level - 1]
            block = acb_mat([[row[i, j] for j in range(2)] for i in range(2)])
            assert block.trace().real > 0
            result[level] = block * (1 / block.trace().real)
        return result

    central = 2_306_888_100_085
    left = projectors(central - 1, (1, 16, 18))
    centre = projectors(central, (7, 10, 16, 18))
    # 单位 Bloch 方向的弦长平方为 4(1 − tr(PP′))。
    distances = {(i, j): 2 * (1 - (p*q).trace().real).sqrt()
                 for i, p in left.items() for j, q in centre.items()}
    pair = min(distances, key=lambda key: distances[key].lower())
    shortest = distances.pop(pair)
    assert all(shortest < value for value in distances.values())
    assert pair == (18, 16) and shortest > arb('.7504203386270') > arb('.748')
    print('投影内积给出的最短层对及弦长', pair, shortest.str(30), flush=True)


def rational_support(forms, spans, box, scale):
    pl, ph, ll, lh, td, rd, loss, sine = [Fraction(int(x), scale) for x in box]
    normals = [[Fraction(int(x), scale) for x in row] for row in forms]
    for angular_forms in spans:
        for (n0, support, a, b, n2), (c, d) in zip(normals, angular_forms, strict=True):
            c, d = Fraction(int(c), scale), Fraction(int(d), scale)
            first = [p*v for p in (pl, ph) for v in (a, b)]
            second = [l*v for l in (ll, lh) for v in (c, d)]
            lo, hi = min(first) + min(second), max(first) + max(second)
            margin = n0 - support/td + lo - max(abs(lo), abs(hi))*loss
            tangent = max(0, n2*rd**2 - max(0, lo, -hi)**2)
            if margin > 0 and margin**2 > tangent*sine**2:
                break
        else:
            return False
    return True


def sample_azimuth():
    fixed = import_module('ch39-integer-azimuth-support')
    fixed.compile_cover(model.ROOT)
    # 两种算法共享相同的向外取整输入。保留两壳的方向表，避免缓存键复用。
    shell_rows = [model.up_rays(offset) for offset in (-1, 0)]
    lower, upper = fixed.lower, fixed.upper
    td = arb('.0002808').union(arb('.000283'))
    r = arb('.9959').union(arb('.9962'))
    parallel = arb('.8122').union(arb('.8132'))
    transverse = arb('.3959').union(arb('.3969'))
    down_radius = arb('.902').union(arb('.906'))
    count = accepted = 0
    for rows in shell_rows:
        for _, _, dist, forms in rows:
            k, delta = 2/dist, 1 - r**2
            assert k*delta < 1
            h = (1-k*delta).sqrt()
            loss = ((k-1)*delta**2 / (r*(r+h)*(1+h)*(1+r))).upper()
            box = (lower(parallel), upper(parallel), lower(transverse), upper(transverse),
                   lower(td), upper(down_radius), upper(loss), upper((2*loss).sqrt()))
            flat, cache = fixed.initialize(forms)
            data = np.array(box, dtype=np.int64)
            pointer = fixed.point
            for depth in (4, 6):
                spans = fixed.section(forms, depth, cache)
                result = fixed.support_library.supports(len(forms), 2**depth, pointer(flat), pointer(spans), pointer(data))
                assert bool(result) == rational_support(flat, spans, box, fixed.M)
                count += 1
                accepted += result
    assert count == 28 and 0 < accepted < count
    print('方位角网格抽查数、两种算法同时判定分离的个数', count, accepted, flush=True)


if __name__ == '__main__':
    replay_with_exact_endpoints()
    chord_distances()
    sample_azimuth()
