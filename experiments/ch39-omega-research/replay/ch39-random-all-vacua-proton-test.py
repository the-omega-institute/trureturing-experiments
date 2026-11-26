"""未知共同真空时，从有限实际编号中随机读一枚的保证均方根。"""
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb


def density_bounds():
    # 正负频率分布在三条整数直线上，精确计算四阶矩中的有序卷积重数。
    support = [(q*a, q*b) for a, b in ((1, 0), (0, 1), (1, -1))
               for q in range(-105, 106) if q]
    counts = Counter((a+x, b+y) for a, b in support for x, y in support)
    assert counts.pop((0, 0)) == len(support) == 630 and max(counts.values()) == 210
    kurtosis = 1 + max(counts.values()) * (1 - Fraction(1, len(support)))
    assert kurtosis == Fraction(632, 3)
    for threshold, claim in ((".0005", ".003874"), (".0008", ".0026876"),
                             (".001", ".00178318"), (".0015", ".0000769166")):
        ratio = (Fraction(threshold)/Fraction(".0016054"))**2
        density = (1-ratio)**2 / (kurtosis-2*ratio+ratio**2)
        assert density > Fraction(claim)
        print("绝对值阈值、实际编号自然密度下界", threshold, claim)


def verify():
    root = Path(__file__).resolve().parent
    detection = run_path(str(root / "ch39-finite-proton-detection.py"))
    experiments = detection["verify"](root)
    density_bounds()
    corners = list(product((0, 1), repeat=3))
    # 表项是非零角点的编号和有理概率的分子，第七角点补足归一。
    shares = ({2: 769736767762, 6: 102783525933991}, {2: 7272023413507})
    for (numbers, filters, upper, deficit), seed in zip(experiments, shares, strict=True):
        mixture = [arb(seed.get(i, 0))/10**15 for i in range(8)]
        mixture[7] = 1 - sum(mixture)
        squares = [[sum(c*v for c, v in zip(row, corner))**2 for corner in corners] for row in filters]
        weights = [sum(lam*v for lam, v in zip(mixture, row)).sqrt() for row in squares]
        assert min(weights) > 0
        lower_norm = sum(weights)
        probabilities = [v/lower_norm for v in weights]
        upper_norm = max(sum(row[j]/p for row, p in zip(squares, probabilities)).sqrt().upper()
                         for j in range(8))
        assert upper_norm-lower_norm < arb("1e-8")
        rms = (1-upper)/upper_norm-deficit
        threshold = ".00016312" if len(numbers) == 218 else ".00036435"
        assert rms > arb(threshold)
        print("实际编号数、最优滤波范数的下上界、单次随机读数的保证均方根",
              len(numbers), lower_norm.str(19), upper_norm.str(19), rms.str(17))


if __name__ == "__main__":
    verify()
