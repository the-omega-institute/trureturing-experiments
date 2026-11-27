"""右外壳保持三循环字符的任意上型输运与任意非负下型权的全域 CKM 障碍。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


TRACE_FLOOR = ".0077"


def transport_bounds(root):
    source = 2_306_888_100_086
    train = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    native = run_path(str(root / "ch39-free-four-site-native-positive-spectrum.py"))["native"]
    frames = run_path(str(root / "ch39-complex-four-site-unit-quark-amplitudes.py"))
    _, labels, _ = train["verify_nodes"](root, ((source, ".0294"),), ".1", (),
                                        ((6, 7, 10, 13), (6, 7, 10, 12)))
    before = native(root, source)[0]
    after, down_after = native(root, source + 1)
    ctx.prec = 384
    rotation = frames["four_basis"]()
    groups = [[], [], []]
    for j, row in enumerate(before):
        rotated = frames["adj"](rotation) * row * rotation
        strengths = [(rotated[0, 0]+rotated[1, 1]).real, rotated[2, 2].real,
                     rotated[3, 3].real]
        group, = [i for i, s in enumerate(strengths) if s > 0]
        groups[group].append(labels[source, -4][j])
    assert [len(g) for g in groups] == [7, 6, 6]
    characters = []
    for group in groups[1:]:
        interval = group[0]
        for value in group[1:]:
            interval = interval.union(value)
        characters.append(interval)
    first, second = sorted(range(2), key=lambda i: characters[i].mid())
    e, f = characters[first], characters[second]

    # 三循环和简单谱确定七条原生平凡线。其余十二条四站投影的平凡块严格为零。
    rows = [row for row in after if (row[0, 0]+row[1, 1]).real > 0]
    assert len(rows) == 7
    coordinates = [((r[0, 0]+r[1, 1]).real, (r[0, 0]-r[1, 1]).real,
                    r[0, 1].real, r[0, 1].imag) for r in rows]
    least = arb(1)
    for weights in permutations(groups[0]):
        trace, difference, real, imag = [
            sum(w*v for w, v in zip(weights, column, strict=True))
            for column in zip(*coordinates, strict=True)
        ]
        minimum = (trace - (difference**2 + 4*(real**2 + imag**2)).sqrt())/2
        least = min(least, (minimum + trace/100).lower())
    assert least > arb(".0000834")
    # 任何给出最轻上型质量的三平面，都要求平凡块最小本征值至多等于 m_u。
    assert e > arb(".00000615")
    assert 100*(least-arb(".00000615")) > arb(TRACE_FLOOR)
    print("字符、七阶置换上的最小支撑量、平凡块迹严格下界",
          e.str(15), f.str(15), least.str(15), (100*(least-arb(".00000615"))).str(15), flush=True)
    down_native = [frames["adj"](rotation)*r*rotation for r in down_after]
    return e, f, down_native, (2+first, 2+second)


def cover(root, e, f, down_native, character_axes):
    pencil = run_path(str(root / "ch39-triple-geometry-pencil.py"))
    lib, geom, mixing, masses, up, down, _ = pencil["setup"](root)
    mu = masses[0]
    gamma = masses[1][2] / mu[2]
    assert mu[0] < e < mu[1] < mu[2] < f
    beta = [(m-e)*(f-m) for m in mu]
    trace_test = [c*(sum(mu)-m-arb(TRACE_FLOOR)) for c, m in zip(beta, mu, strict=True)]
    det_test = [c*mu[(i+1)%3]*mu[(i+2)%3] for i, c in enumerate(beta)]

    # 展示一个两荷四站扩张均正定的代数点，以及此点与原生上型迹界的距离。
    t = arb(".1")
    y = lib["lines"](up, down, t, t, (2,), geom)[2]
    s, up_prob, _ = lib["moments"](mixing*y, mu)
    sd, down_prob, _ = lib["moments"](y, masses[1])
    ed, fd = [sd+gamma*(s-v)/t for v in (e, f)]
    assert min(ed, fd) > 0
    for q, (m, p, left, right) in enumerate(((mu, up_prob, e, f),
                                              (masses[1], down_prob, ed, fd))):
        weights = [pi*(mi-left)*(right-mi) for pi, mi in zip(p, m, strict=True)]
        total = sum(weights)
        assert total > 0
        trace = sum(w*(sum(m)-mi) for w, mi in zip(weights, m, strict=True))/total
        determinant = sum(w*m[(i+1)%3]*m[(i+2)%3] for i, w in enumerate(weights))/total
        gap = (trace**2-4*determinant).sqrt()
        small, large = (trace-gap)/2, (trace+gap)/2
        assert 0 < small < large
        print("参数 0.1 的部门、两字符、平凡块本征值", q, left.str(12), right.str(12),
              small.str(12), large.str(12), flush=True)
        if q == 0:
            assert trace < arb(".003117")
        else:
            # 在下型原生坐标中选第一条平凡射线作较大本征线。
            left_axis, right_axis = character_axes
            trivial = [j for j,r in enumerate(down_native) if (r[0,0]+r[1,1]).real > 0]
            chosen = trivial[0]
            weights = [small if j in trivial else
                       left if r[left_axis,left_axis].real > 0 else right
                       for j,r in enumerate(down_native)]
            row = down_native[chosen]
            weights[chosen] += (large-small)/(row[0,0]+row[1,1]).real
            assert len(trivial) == 7 and min(weights) > 0
            raw = sum((w*r for w,r in zip(weights, down_native, strict=True)), geom["diag"]([0]*4))
            assert raw[left_axis,left_axis].real.overlaps(left)
            assert raw[right_axis,right_axis].real.overlaps(right)
            assert (raw[0,0]+raw[1,1]).real.overlaps(trace)
            assert (raw[0,0]*raw[1,1]-raw[0,1]*raw[1,0]).real.overlaps(determinant)
            print("下型十九权最小、最大", min(weights).str(12), max(weights).str(12), flush=True)

    def excluded(y, chart, t, sign):
        s, probability, _ = lib["moments"](mixing*y, mu)

        def weighted(coefficients):
            return sum(p*c for p, c in zip(probability, coefficients, strict=True))

        if s < e or s > f or weighted(beta) < 0:
            return True
        if weighted(trace_test) < 0 or weighted(det_test) < 0:
            return True
        sd = lib["moments"](y, masses[1])[0]
        if chart == 1:
            return any(sd + gamma*t*(s-v) < 0 for v in (e, f))
        # 这一张图中先乘 t 的绝对值，零附近也保留所有非零 t 的可行点。
        return any(sign*(t*sd + gamma*(s-v)) < 0 for v in (e, f))

    for chart, (base, other) in enumerate(((up, down), (down, up))):
        for branch in range(3):
            stack = [(arb(-1), arb(0), -1, 0), (arb(0), arb(1), 1, 0)]
            leaves = deepest = 0
            while stack:
                lo, hi, sign, depth = stack.pop()
                lines = lib["lines"](base, other, lo, hi, (branch,), geom)
                if branch in lines and excluded(lines[branch], chart, lo.union(hi), sign):
                    leaves += 1
                    deepest = max(deepest, depth)
                    continue
                assert depth < 25, (chart, branch, lo, hi)
                mid = (lo+hi)/2
                stack.extend(((lo, mid, sign, depth+1), (mid, hi, sign, depth+1)))
            print("图、支、覆盖叶区间、最大深度", chart, branch, leaves, deepest, flush=True)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    cover(root, *transport_bounds(root))
