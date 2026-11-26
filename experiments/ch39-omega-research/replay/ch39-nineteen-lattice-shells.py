"""完整非负八分体整数球壳的十九点分类与两种相似类的距离证书。"""

from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import isqrt


def shell(norm):
    points = []
    for x in range(isqrt(norm) + 1):
        for y in range(isqrt(norm - x * x) + 1):
            remaining = norm - x * x - y * y
            z = isqrt(remaining)
            if z * z == remaining:
                points.append((x, y, z))
    return points


def odd_local_factor(p, k):
    # χ_{−3}(p) 在 3 处为零，在其余奇素数上由模 3 余数决定。
    character = (0, 1, -1)[p % 3]
    return 1 + (p - character) * sum(p**j for j in range(k))


def check_count():
    # 局部因子至少为 p，故只剩以下奇素数。对 p=3 检查至首次超过 19。
    primes = [p for p in range(3, 20, 2) if all(p % d for d in range(2, isqrt(p) + 1))]
    smallest = sorted(odd_local_factor(p, 1) for p in primes)
    assert smallest[0] * smallest[1] == 28 > 19
    candidates = {1: 1}
    for p in primes:
        k = 1
        while odd_local_factor(p, k) <= 19:
            candidates[p**k] = odd_local_factor(p, k)
            k += 1
    assert candidates == {1: 1, 3: 4, 5: 7, 7: 7, 9: 13, 11: 13, 13: 13, 17: 19, 19: 19}
    for t, count in sorted(candidates.items()):
        points = shell(3 * t**2)
        distances = {square_distance(p, q) for p, q in combinations(points, 2)}
        assert len(points) == count
        assert all(d % 4 == 0 for d in distances)
        numerators = sorted(d for d in distances if d % 181 == 0)
        assert numerators == ([724] if t == 17 else [])
        if t == 17:
            assert 32 in distances
        print("小壳的奇因子、点数、最大平方弦长、181 的倍数", t, count,
              max(distances, default=0), numerators)


def square_distance(first, second):
    return sum((a - b)**2 for a, b in zip(first, second))


def volume6(center, anchors):
    u, v, w = [tuple(a - b for a, b in zip(point, center)) for point in anchors]
    return sum(u[i] * (v[(i + 1) % 3] * w[(i + 2) % 3]
                       - v[(i + 2) % 3] * w[(i + 1) % 3]) for i in range(3))


def compare(t):
    points = shell(3 * t**2)
    assert len(points) == len(set(points)) == 19
    center = (t, t, t)
    orbits = Counter(tuple(sorted(p)) for p in points)
    radius = {p: square_distance(center, p) for p in points}
    eccentricity = {p: max(square_distance(p, q) for q in points) for p in points}
    rows = sorted((p, size, radius[p], eccentricity[p]) for p, size in orbits.items())
    assert sorted(orbits.values()) == [1, 3, 3, 6, 6]
    assert len({eccentricity[p] for p in orbits}) == len(orbits)
    assert [p for p in points if eccentricity[p] == min(eccentricity.values())] == [center]
    print("t、排序轨道代表、轨道大小、到中心平方距、最大平方距", t, rows)

    # 每个轨道均由最大平方距辨认。一个三点轨道与中心仿射无关。
    anchor = (13, 13, 23) if t == 17 else (11, 11, 29)
    tetrahedron = sorted(set(permutations(anchor)))
    assert set(tetrahedron) == {p for p in points if eccentricity[p] == eccentricity[anchor]}
    assert abs(volume6(center, tetrahedron)) == (200 if t == 17 else 1944)
    # 保持四点的等距由三点置换决定。全部六种坐标置换在整壳上实现。
    for permutation in permutations(range(3)):
        action = {tuple(p[i] for i in permutation) for p in points}
        assert action == set(points)

    near_radius = min(d for d in radius.values() if d > 0)
    far_radius = max(radius.values())
    near = [p for p in points if radius[p] == near_radius]
    far = [p for p in points if radius[p] == far_radius]
    far_pairs = list(combinations(far, 2))
    base = min(square_distance(p, q) for p, q in far_pairs)
    bases = [(p, q) for p, q in far_pairs if square_distance(p, q) == base]
    assert len(bases) == 3 and {p for pair in bases for p in pair} == set(far)
    expected_sides = ({(244, 244): 1, (724, 804): 2} if t == 17 else
                      {(216, 296): 2, (536, 776): 2, (776, 936): 2})
    for p, q in bases:
        sides = Counter(tuple(sorted((square_distance(p, z), square_distance(q, z)))) for z in near)
        assert sides == expected_sides
    assert (len(near), len(far), near_radius, far_radius, base) == (
        (3, 6, 68, 544, 32) if t == 17 else (6, 6, 76, 532, 200))
    print("最近和最远点数、平方半径比、底边平方长、各侧边对",
          len(near), len(far), Fraction(near_radius, far_radius), base, expected_sides)

    distances = Counter(square_distance(p, q) for p, q in combinations(points, 2))
    target = Fraction(181, 8)
    assert distances.total() == 171
    pairs = [(s, b) for s, b in product(distances, repeat=2) if Fraction(s, b) == target]
    assert pairs == ([(724, 32)] if t == 17 else [])
    ratios = {Fraction(a, b) for a, b in combinations(sorted(distances, reverse=True), 2)}
    assert (len(distances), len(ratios)) == ((30, 412) if t == 17 else (27, 316))
    gap = min(abs(r - target) for r in ratios)
    closest = sorted(r for r in ratios if abs(r - target) == gap)
    assert (gap, closest) == ((0, [target]) if t == 17 else (Fraction(37, 56), [Fraction(163, 7)]))
    print("全部平方弦长", sorted(distances))
    print("平方弦长比等于 181/8 的有序长度对、最近比值、绝对差", pairs, closest, gap)
    return distances


def compare_original_interval(distances):
    # 只与第一族的菜单和数值区间作算术比较。第二族自身的曲率界依赖其层谱。
    lower, upper = Fraction("31.9830"), Fraction("32.0170")
    records = [(Fraction(724 * b, a)**2 / k, k, a, b)
               for k, a, b in product((16, 32, 48, 96), distances, distances)]
    below = max(record for record in records if record[0] < lower)
    above = min(record for record in records if record[0] > upper)
    passing = [record for record in records if lower <= record[0] <= upper]
    assert len(records) == 2916 and passing == []
    assert below == (Fraction(1605289, 50625), 16, 1800, 56)
    assert above == (Fraction(11826721, 367500), 48, 1400, 76)
    print("t=19 套用第一族数值菜单的候选数、落入第一族区间的项数、两侧最近项",
          len(records), len(passing), below, above)


if __name__ == "__main__":
    check_count()
    compare(17)
    compare_original_interval(compare(19))
