"""在 11431 的孤立中点窗口中核查最优正权的碰撞占比。"""
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    ctx.prec = 192
    start = 11422
    heights = [z.imag for z in acb.zeta_zeros(start, 19)]
    pairs = list(combinations_with_replacement(range(19), 2))
    sums = [heights[i] + heights[j] for i, j in pairs]
    gaps = [(abs(sums[a] - sums[b]), pairs[a], pairs[b])
            for a, b in combinations(range(len(pairs)), 2)]
    previous = run_path(str(Path(__file__).with_name("ch39-zero-midpoint-fejer.py")))
    first, second = previous["separated_first_two"](gaps)
    assert first[0] > 0 and first[1:] == ((5, 14), (9, 9))
    delta, radius = first[0], second[0]
    ramp = 1 - delta / radius

    # 权重按左端、右端、中点的次序为 1、1、6^(1/4)，再归一。
    relative_center = arb(6).sqrt().sqrt()
    end_weight = 1 / (2 + relative_center)
    center_weight = relative_center * end_weight
    weights = {5: end_weight, 14: end_weight, 9: center_weight}
    occupied = [((i, j), (1 if i == j else 2) * weights[i] * weights[j])
                for i, j in pairs if i in weights and j in weights]
    floor = sum(mass**2 for _, mass in occupied)
    excess = arb(0)
    for ((i, j), a), ((k, l), b) in combinations(occupied, 2):
        gap = abs(heights[i] + heights[j] - heights[k] - heights[l])
        excess += 2 * a * b * (1 - gap / radius).max(0)
    actual_share = excess / (floor + excess)
    sharp_share = 4 * ramp / (8 + 2 * arb(6).sqrt() + 4 * ramp)
    hard_share = 4 / (12 + 2 * arb(6).sqrt())
    assert actual_share.overlaps(sharp_share)
    assert arb(".23657016") < actual_share < arb(".23657018")
    assert arb(".23670068") < hard_share < arb(".23670069")
    print("唯一额外配对、端点与中点的最优权重、最大三角核占比、硬带占比",
          (start + 5, start + 14, start + 9), end_weight.str(12), center_weight.str(12),
          sharp_share.str(14), hard_share.str(14))


if __name__ == "__main__":
    verify()
