"""逐零点检验十个共同约束，并用支撑列与约化成本认证 CP 读数的端点。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

previous = run_path(str(Path(__file__).with_name("ch39-eleven-neutral-cp-range.py")))
joint = previous["joint"]
extra = (142270888235, 142269938072, 142271215525, 142270381883,
         295874150547, 295872524995, 295874398386, 295872666901)
first_list = (*previous["indices"], *extra)
high_indices = (143921294494500, 143921293719910, 143921298450585, 143921293614265)
indices = (*first_list, *high_indices)
lower_support = (8359, 143921293989490, 295874398386, 142270715263, 13219432687559,
                 295668865161, 142270888235, 142269606336, 142268163595, 295873886415)
upper_support = (8359, 143921293989490, 295872666901, 142271930613, 13219432687559,
                 295668865161, 295873014395, 142269606336, 142268163595, 295872524995)

high_upper_support = (8359, 295874398386, 142271930613, 142270381883, 13219432687559,
                      295668865161, 295872666901, 295872524995, 142268163595, 143921293614265)


def endpoint(indices, rows, target, columns, support, direction, invariants, ref, current):
    members = [indices.index(n) for n in support]
    matrix = arb_mat(list(zip(*(columns[i] for i in members), strict=True)))
    weights = matrix.solve(arb_mat([[value] for value in target])).entries()
    assert all(w > 0 for w in weights)
    flavor = [rows[i][2]*w for i, w in zip(members, weights)]
    total = sum(flavor)
    joint["replay"]["neutral"]["verify_maximum"](
        [rows[i][1] for i in members], [v/total for v in flavor], invariants, ref)
    single = [1-r[4]/current for r in rows]
    loss = sum(v*single[i] for i, v in zip(members, flavor))/total

    # 将其余各列在该端点的十列基上精确展开，检验有向目标的约化成本。
    objective = [direction*r[2]*(h-loss) for r, h in zip(rows, single)]
    for j, column in enumerate(columns):
        if j in members:
            continue
        coordinates = matrix.solve(arb_mat([[value] for value in column])).entries()
        reduced = objective[j] - sum(objective[i]*v for i, v in zip(members, coordinates))
        assert reduced < 0
    return loss


def problem(rows, times, probabilities):
    columns = [[*r[0], r[2]*(r[3][7]-probabilities[7]-arb("2e-12")),
                joint["model"]["time_kernel"](r[1], 6)] for r in rows]
    return columns, [*times[:6], 0, 0, 0, times[6]+arb("1.360900010e-8")]


def verify():
    rows, times, invariants, ref, probabilities, current = joint["columns_at"](2_903_600_000, indices)
    ctx.prec = 640
    columns, target = problem(rows, times, probabilities)
    lo = endpoint(indices, rows, target, columns, lower_support, -1, invariants, ref, current)
    old_size = len(first_list)
    old_hi = endpoint(first_list, rows[:old_size], target, columns[:old_size], upper_support,
                      1, invariants, ref, current)
    hi = endpoint(indices, rows, target, columns, high_upper_support, 1, invariants, ref, current)
    assert arb(".05570022060867") < lo < arb(".05570022060869")
    assert arb(".07015185002496") < old_hi < arb(".07015185002499")
    assert arb(".093435233874") < hi < arb(".093435233877")
    print("十九枚、二十三枚共用的下端和各自的上端", lo.str(17), old_hi.str(17), hi.str(15))


if __name__ == "__main__":
    verify()
