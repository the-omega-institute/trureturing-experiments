"""三窗、五个算术时间及逐背景夸克概率的联合读出。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

root = Path(__file__).parent
replay = run_path(str(root / "ch39-seven-actual-neutral-times.py"))
cp = run_path(str(root / "ch39-cp-ensemble-test.py"))
flavor = cp["flavor"]
model = replay["neutral"]["model"]
eight_indices = (8359, 143921293989490, 295873368963, 142270715263,
                 13219432687559, 295668865161, 295873014395, 142269606336)


def columns_at(scale, indices):
    ctx.prec, ctx.cap = 768, 3
    times = model["shared"]["clock"]["clock"]["times"]()
    invariants, ref, poles, slopes = model["shared"]["flavor_data"]()
    matrix = flavor["reference_matrix"]()
    phase, _ = flavor["singlet"](11661)
    prob0, _, cp0 = cp["matrix_readings"](matrix, phase)
    records = []
    for height, x, eta in replay["shell_nodes"](scale, indices, height_bits=176):
        clocks = flavor["geom"]["phases"](height)
        phase, _ = flavor["singlet_at"](clocks)
        probabilities, _, cp_value = cp["matrix_readings"](matrix, phase)
        times_at = [model["time_kernel"](x, j) for j in range(1, 6)]
        column = [1, *times_at, eta*model["gradient"](x, poles, slopes), eta*(probabilities[6]-prob0[6])]
        records.append((column, x, eta, probabilities, cp_value))
    return records, times, invariants, ref, prob0, cp0


def solve(records, times):
    matrix = arb_mat(list(zip(*(row[0] for row in records), strict=True)))
    target = arb_mat([[v] for v in [*times[:6], 0, 0]])
    weights = matrix.solve(target)
    assert all(v.contains(0) for v in (matrix*weights-target).entries())
    return weights.entries()


def compare(records, weights, times, invariants, ref, prob0, cp0):
    raw = [w*row[2] for w,row in zip(weights, records, strict=True)]
    total = sum(raw)
    flavor_weights = [v/total for v in raw]
    nodes = [row[1] for row in records]
    replay["neutral"]["verify_maximum"](nodes, flavor_weights, invariants, ref)
    probs = [sum(w*row[3][i] for w,row in zip(flavor_weights,records,strict=True)) for i in range(9)]
    cp_mean = sum(w*row[4] for w,row in zip(flavor_weights, records, strict=True))
    sixth = sum(w*model["time_kernel"](x, 6) for w,x in zip(weights,nodes,strict=True))
    print("算术概率",*[w.str(10) for w in weights])
    print("平均 |Vub|² 的相对误差", (probs[6]/prob0[6]-1).str(12))
    print("平均 |Vcb|² 的绝对偏差", (probs[7]-prob0[7]).str(12))
    print("平均 J 的相对亏额", (1-cp_mean/cp0).str(12))
    print("第六项算术时间偏高", (sixth-times[6]).str(12))
    return probs[6]-prob0[6], probs[7]-prob0[7], 1-cp_mean/cp0


def verify():
    records, times, invariants, ref, probabilities, j0 = columns_at(2_903_600_000, eight_indices)
    weights = solve(records, times)
    assert all(w > 0 for w in weights)
    ub, cb, cp_loss = compare(records, weights, times, invariants, ref, probabilities, j0)
    assert ub.contains(0) and arb("1.90254e-12") < cb < arb("1.90255e-12")
    assert arb(".06083301837") < cp_loss < arb(".06083301839")


if __name__ == "__main__":
    verify()
