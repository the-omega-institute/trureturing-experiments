"""在十九站整数壳上以配对弦距比的有理性筛出三条极小弧。"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import isqrt
from pathlib import Path
from runpy import run_path


def inspect(sites, radius):
    counts, accepted = Counter(), []
    for (i,a), (j,b) in combinations(enumerate(sites), 2):
        d = sum((x-y)**2 for x,y in zip(a,b,strict=True))
        counts[d] += 1
        two_kappa = Fraction(4*radius - d, 2*radius)
        assert two_kappa > 0
        if (isqrt(two_kappa.numerator)**2 == two_kappa.numerator
                and isqrt(two_kappa.denominator)**2 == two_kappa.denominator):
            accepted.append((i,j,d))
    return counts, accepted


def quadratic_square_class(radius, distance):
    """正有理数 2κ 的平方类用唯一的无平方因子正整数表示。"""
    value = Fraction(4*radius-distance, 2*radius)
    n = value.numerator * value.denominator
    factor = 2
    while factor*factor <= n:
        square = factor*factor
        while n % square == 0:
            n //= square
        factor += 1
    return n


def verify(root):
    sites = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))["geo"]["sites"]
    assert len(sites) == 19 and all(sum(x*x for x in p) == 867 for p in sites)
    counts, accepted = inspect(sites, 867)
    assert sum(counts.values()) == 171 and len(counts) == 30
    assert accepted == [(6,11,12),(7,12,12),(8,18,12)]
    assert min(counts) == 12 and counts[12] == 3
    print("无序站对数、不同平方弦长数、有理弦距比候选",sum(counts.values()),len(counts),accepted)

    for radius, expected_count in ((17, 9), (867, 19), (1083, 19)):
        whole = [p for p in product(range(isqrt(radius)+1), repeat=3)
                 if sum(x*x for x in p) == radius]
        assert len(whole) == expected_count
        lengths, selected = inspect(whole, radius)
        if radius == 867:
            assert set(whole) == set(sites) and len(selected) == 3
        else:
            assert not selected
        print("完整非负整数壳的平方半径、点数、通过检验的站对数",radius,len(whole),len(selected))
        if radius in (867, 1083):
            classes = {d: quadratic_square_class(radius, d) for d in lengths}
            assert len(set(classes.values())) == len(lengths)
            visible = {d: n for d, n in lengths.items() if classes[d] in (1, 5154)}
            assert visible == ({12: 3, 32: 3} if radius == 867 else {})
            print("在有理数域添根号5154后的平方弦长及条数",radius,visible)
    assert Fraction(4*867-32, 2*867) == Fraction(5154, 51**2)



if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
