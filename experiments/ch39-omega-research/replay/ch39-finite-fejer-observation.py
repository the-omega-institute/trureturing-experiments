"""十九枚实际零点的有限时间四阶 Fejér 积分首次超过占据对角。"""
from itertools import combinations, combinations_with_replacement

from flint import acb, arb, ctx


def finite_transform(gap, resolution, radius):
    """偶核在 [-R,R] 上的余弦积分，由分部积分精确给出。"""
    def primitive(x):
        z = radius * x
        return x * z.si() - 2 * (z / 2).sin() ** 2 / radius

    return (primitive(gap + resolution) + primitive(gap - resolution)
            - 2 * primitive(gap)) / (arb.pi() * resolution)


def verify():
    ctx.prec = 160
    eta = arb(1) / 1000
    heights = [z.imag for z in acb.zeta_zeros(11422, 19)]
    pairs = list(combinations_with_replacement(range(len(heights)), 2))
    multiplicities = [1 if i == j else 2 for i, j in pairs]
    pair_sums = [heights[i] + heights[j] for i, j in pairs]
    diagonal = sum(m * m for m in multiplicities)
    assert diagonal == 703
    differences = []
    isolated = []
    for i, j in combinations(range(len(pairs)), 2):
        gap = abs(pair_sums[i] - pair_sums[j])
        weight = 2 * multiplicities[i] * multiplicities[j]
        differences.append((gap, weight))
        if gap < eta:
            isolated.append((pairs[i], pairs[j], gap, weight))
        else:
            assert gap > eta
    (left, right, gap, weight), = isolated
    assert (left, right, weight) == ((5, 14), (9, 9), 4)
    excess = weight * (1 - gap / eta)
    sufficient = 64 * (diagonal + excess) / (eta * excess)
    assert 11330942 < sufficient < 11330943
    print("通用截短界保证超过对角的最小整数半径", 11330943)

    for radius, lower, upper in [(111520, "-0.0000025", "-0.0000024"),
                                 (111521, "0.000038", "0.000039")]:
        total = diagonal * finite_transform(arb(0), eta, radius)
        for distance, factor in differences:
            total += factor * finite_transform(distance, eta, radius)
        margin = total - diagonal
        assert arb(lower) < margin < arb(upper)
        print("观察半径与乘以 19^4 的超额", radius, margin.str(14))


if __name__ == "__main__":
    verify()
