"""用零点区间比较三千个中点，并分离十九点 Fejér 矩中的占据对角。"""
from itertools import combinations, combinations_with_replacement

from flint import acb, arb, ctx


def separated_first_two(rows):
    # 浮点仅提出候选；严格区间比较决定前两名。
    candidates = sorted(range(len(rows)), key=lambda k: float(rows[k][0]))
    first, second = candidates[:2]
    assert rows[first][0] < rows[second][0]
    assert all(rows[second][0] < rows[k][0] for k in candidates[2:])
    return rows[first], rows[second]


def midpoint_scan(heights, start):
    rows = []
    for n in range(10000, 13000):
        k = n-start
        twice = 2*heights[k]
        options = [(abs(heights[k-i] + heights[k+j] - twice), i, j)
                   for i in range(1, 10) for j in range(1, 10)]
        # 每个中心的 81 个候选均与其报出的最小值严格分离。
        pos = min(range(len(options)), key=lambda idx: float(options[idx][0]))
        error, i, j = options[pos]
        assert all(error < item[0] for idx, item in enumerate(options) if idx != pos)
        spacing = (heights[k+9] - heights[k-9])/18
        rows.append((error, n, i, j, error/spacing))
    first, second = separated_first_two(rows)
    assert first[1:4] == (11431, 4, 5) and second[1] == 11628
    assert arb("35.59") < second[0]/first[0] < arb("35.60")
    normalized = [(r[4], r[1]) for r in rows]
    assert separated_first_two(normalized)[0][1] == 11431
    for offset, expected in ((0, 10074), (1000, 11431), (2000, 12783)):
        assert separated_first_two(rows[offset:offset+1000])[0][1] == expected
    print("三千中心前两名、两倍中点偏差及其比值",
          first[1:4], first[0].str(14), second[1:4], second[0].str(14),
          (second[0]/first[0]).str(12))


def compare_window(heights):
    points = list(combinations_with_replacement(range(19), 2))
    pairs = [(heights[i]+heights[j], 1 if i == j else 2) for i, j in points]
    diagonal = sum(m*m for _, m in pairs)
    assert diagonal == 703 and sum(m for _, m in pairs) == 19**2
    differences = [(abs(pairs[i][0]-pairs[j][0]), 2*pairs[i][1]*pairs[j][1],
                    points[i], points[j]) for i, j in combinations(range(len(pairs)), 2)]
    nearest, runner_up = separated_first_two(differences)
    assert nearest[2:] == ((5,14), (9,9)) and nearest[1] == 4
    assert arb("0.0000016812114283") < nearest[0] < arb("0.0000016812114285")
    assert arb("0.0023277388") < runner_up[0] < arb("0.0023277390")
    radius = arb(1)/100
    inside = [d for d in differences if d[0] < radius]
    assert all(d[0] < radius or d[0] > radius for d in differences)
    assert len(inside) == 10 and sum(d[1] for d in inside) == 72
    excess = sum(m*(1-gap/radius) for gap, m, _, _ in inside)
    moment = (diagonal + excess)/19**4
    share = diagonal/(diagonal + excess)
    assert arb(".00570094240512") < moment < arb(".00570094240513")
    assert arb(".94622467201") < share < arb(".94622467203")
    print("最近两种不同配对的距离、Fejér 四阶矩、占据对角占比",
          nearest[0].str(14), runner_up[0].str(14), moment.str(14), share.str(14))


def dyadic_coefficients(heights):
    # 临界线上的 rho，使冻结不等式的右移系数成为以下实函数。
    log_two = arb(2).log()
    coeff = [(3-2*arb(2).sqrt()*(t*log_two).cos())/(2*(t*t+arb(1)/4))
             for t in heights]
    order = sorted(range(len(coeff)), key=lambda k: -float(coeff[k]))
    assert all(coeff[a] > coeff[b] for a, b in zip(order, order[1:]))
    assert order[0] == 0 and order[3] == 9
    assert arb("2.20461e15") < 1/coeff[9]**2 < arb("2.20462e15")
    print("倍频莫比乌斯壳下界系数、11431 名次、常数一交点",
          coeff[9].str(14), order.index(9)+1, (1/coeff[9]**2).str(14))


if __name__ == "__main__":
    ctx.prec = 160
    start = 9991
    heights = [z.imag for z in acb.zeta_zeros(start, 3018)]
    assert len(heights) == 3018 and all(x < y for x, y in zip(heights, heights[1:]))
    window = heights[11422-start:11441-start]
    compare_window(window)
    dyadic_coefficients(window)
    midpoint_scan(heights, start)
