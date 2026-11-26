"""固定五项算术时间、一个或两个味概率时，十一枚实际零点的 CP 可行区间。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

joint = run_path(str(Path(__file__).with_name("ch39-five-neutral-quark-probability.py")))
indices = (*joint["eight_indices"], 142268163595, 295873886415, 142271930613)


def cp_loss(rows, weights, reference):
    flavor = [w*row[2] for w, row in zip(weights, rows, strict=True)]
    return sum(v*(1-row[4]/reference) for v, row in zip(flavor, rows, strict=True))/sum(flavor)


def feasible_vertices(rows, times):
    for support in combinations(range(len(rows)), 8):
        # 该尺度的五项时间和驻点强迫第 8359 枚带严格正权。
        if support[0] != 0:
            continue
        subset = [rows[i] for i in support]
        weights = joint["solve"](subset, times)
        if any(v < 0 for v in weights):
            continue
        assert all(v > 0 for v in weights)
        yield support, subset, weights


def five_time_range(rows, times, invariants, phase, reference):
    losses = []
    for _, subset, weights in feasible_vertices(rows, times):
        flavor = [w*r[2] for w, r in zip(weights, subset, strict=True)]
        joint["replay"]["neutral"]["verify_maximum"](
            [r[1] for r in subset], [v/sum(flavor) for v in flavor], invariants, phase)
        losses.append(cp_loss(subset, weights, reference))
    assert len(losses) == 10
    ordered = sorted(losses, key=float)
    assert all(a < b for a, b in zip(ordered, ordered[1:]))
    assert arb(".05828206184780") < ordered[0] < arb(".05828206184783")
    assert arb(".07601188899775") < ordered[-1] < arb(".07601188899779")
    print("十一枚的顶点数及平均 J 的相对亏额两端", len(losses),
          ordered[0].str(16), ordered[-1].str(16))


def pinned_range(rows, times, charm_reference, reference):
    charm_shift = arb("2e-12")
    sixth_shift = arb("1.360900010e-8")
    columns = [[*column, eta*(probs[7]-charm_reference-charm_shift),
                joint["model"]["time_kernel"](x, 6)] for column, x, eta, probs, _ in rows]
    target = arb_mat([[v] for v in [*times[:6], 0, 0, 0, times[6]+sixth_shift]])
    # 以第十一枚的算术权为 t，前十个权由同一个可逆矩阵决定。
    basis = arb_mat(list(zip(*columns[:10], strict=True)))
    start = basis.solve(target).entries()+[arb(0)]
    velocity = basis.solve(arb_mat([[-v] for v in columns[10]])).entries()+[arb(1)]
    assert all(v > 0 for v in start[:10])
    exit_index = 3
    assert velocity[exit_index] < 0
    end_time = -start[exit_index]/velocity[exit_index]
    end = [a+end_time*v for a, v in zip(start, velocity, strict=True)]
    assert 0 < end_time < 1 and end[exit_index].contains(0)
    assert all(v > 0 for i, v in enumerate(end) if i != exit_index)
    lo, hi = [cp_loss(rows, weights, reference) for weights in (start, end)]
    assert arb(".06169080517549") < lo < arb(".06169080517551")
    assert arb(".06413182131993") < hi < arb(".06413182131996")
    print("另固定第六项和平均 |Vcb|² 后，线段长度参数和 J 亏额两端",
          end_time.str(16), lo.str(16), hi.str(16))


def verify():
    rows, times, invariants, phase, probabilities, reference = joint["columns_at"](2_903_600_000, indices)
    ctx.prec = 640
    five_time_range(rows, times, invariants, phase, reference)
    pinned_range(rows, times, probabilities[7], reference)


if __name__ == "__main__":
    verify()
