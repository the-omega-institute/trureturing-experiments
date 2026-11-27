"""在十九站整数壳上以配对弦距比的有理性筛出三条极小弧。"""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import isqrt
from pathlib import Path
from runpy import run_path


def verify(root):
    sites = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))["geo"]["sites"]
    radius = 867
    assert len(sites) == 19 and all(sum(x*x for x in p) == radius for p in sites)
    counts, accepted = Counter(), []
    for (i,a), (j,b) in combinations(enumerate(sites), 2):
        d = sum((x-y)**2 for x,y in zip(a,b,strict=True))
        counts[d] += 1
        two_kappa = Fraction(4*radius - d, 2*radius)
        assert two_kappa > 0
        if (isqrt(two_kappa.numerator)**2 == two_kappa.numerator
                and isqrt(two_kappa.denominator)**2 == two_kappa.denominator):
            accepted.append((i,j,d))
    assert sum(counts.values()) == 171 and len(counts) == 30
    assert accepted == [(6,11,12),(7,12,12),(8,18,12)]
    assert min(counts) == 12 and counts[12] == 3
    print("无序站对数、不同平方弦长数、有理弦距比候选",sum(counts.values()),len(counts),accepted)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
