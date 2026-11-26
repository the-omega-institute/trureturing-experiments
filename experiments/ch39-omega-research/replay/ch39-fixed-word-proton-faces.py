"""用同一固定真空权比较原生面、按类加权与逐面加权的质子树级因子。"""
from collections import Counter
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def signals(root):
    def load(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    faces = load("mincut-real-triangle-tomography")
    shell = load("three-ray-chiral-boundary")
    lib = load("mincut-global-winding")
    fixed = load("height-independent-flavor-vacuum")
    model = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    ctx.prec = 256
    x = arb("-.01046")
    y = x - model["reference_angles"](flavor)[0]
    weights = fixed["minimum_order_weights"](lib, lib["terms"](shell), x, y)
    height = acb.zeta_zero(2_306_888_100_085).imag
    delta = arb(23)/25 * arb(2)**-20 / (height**2 + shell["soft"]["point"]["study"]["offset"])
    phases = shell["geo"]["phases"](height)
    actual = (phases[724].arg(), phases[804].arg())
    pole, tuned = shell["soft"]["pole"], shell["soft"]["star"]["orbit"]
    types = Counter()
    for triangle, (clock, phase) in faces["enumerate_faces"]().items():
        if clock != (0, 0):
            soft = int(pole in triangle and any(k in tuned for k in triangle))
            types[clock, phase, soft] += 1
    assert len(types) == 32 and set(types.values()) == {3}
    rows = []
    for (a, b), pair, soft in sorted(types):
        z = faces["unit"](pair)
        def read(xx, yy):
            loop = z * acb(0, a*xx+b*yy).exp()
            return sum(w * (loop**q).real for q, w in weights.items())
        ref = read(x, y)
        num = read(*actual) * (1 - delta*soft)
        rows.append((ref, ref*soft, (num-ref)*10**6))
    return sorted(types), rows, delta


def minimum_variation(columns, forced, basis, positive):
    # 三列约束分别给出总质量及两个质子因子条件。
    n = len(columns)
    uniform = arb(1) / n
    matrix = arb_mat([[columns[j][i] for j in basis] for i in range(3)])
    residual = arb_mat([[int(i == 0) - sum(col[i] * uniform for j, col in enumerate(columns)
                                            if j not in forced)] for i in range(3)])
    change = matrix.solve(residual).entries()
    dual = matrix.transpose().solve(arb_mat([[int(j in positive)] for j in basis])).entries()
    for j, x in zip(basis, change):
        if j in positive:
            assert x > 0
        else:
            assert -uniform < x < 0
    for j, col in enumerate(columns):
        if j in forced or j in basis:
            continue
        score = sum(y * a for y, a in zip(dual, col))
        assert 0 < score < 1
    tv = sum(x for j, x in zip(basis, change) if j in positive)
    assert (tv - sum(y * b for y, b in zip(dual, residual.entries()))).contains(0)
    print("清零轨道、调整轨道、概率增量、最小总变差",
          forced, basis, [v.str(16) for v in change], tv.str(16))
    return tv


def verify():
    types, rows, delta = signals(Path(__file__).resolve().parent)
    cap = arb(23) / 25 * arb(2)**-20 / 14**2
    for clock in ((1, 0), (0, 1), (1, -1)):
        subset = [row for key, row in zip(types, rows) if key[0] == clock]
        at_reference = sum(v - arb(0).union(cap) * soft for v, soft, diff in subset) / len(subset)
        at_number = sum(v + diff/10**6 for v, soft, diff in subset) / len(subset)
        assert at_reference > arb(".03143") and at_number > arb(".03143")
        print("类别、目标角全高度及实际零点的因子", clock, at_reference.str(13), at_number.str(13))
    finite = [(arb(1), ref-delta*soft, diff+10**6*delta*soft) for ref, soft, diff in rows]
    t1 = minimum_variation(finite, (), (2, 4, 23), (2,))
    assert arb(".03544040262") < t1 < arb(".03544040264")
    forced = tuple(j for j, key in enumerate(types) if key[2])
    assert forced == (2, 21) and all(rows[j][1] < 0 for j in forced)
    all_heights = [(arb(1), ref, diff) for ref, soft, diff in rows]
    t2 = minimum_variation(all_heights, forced, (10, 27, 23), (10, 27))
    assert arb(".08581786483") < t2 < arb(".08581786485")
    for i in (2, 4, 10, 21, 23, 27):
        print("调整轨道的索引及几何标签", i, types[i])


if __name__ == "__main__":
    verify()
