"""有限偶阶计数的完整角亏误差、精确平直驻点与交叉偏振色散。"""
from fractions import Fraction
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
import sympy as sp


def mean_defect(n, r):
    """返回 1-E[k/N]/r。交替余项和正项尾界避开两个指数小量的消差。"""
    y = 2*n*r
    last = y**(2*n)*(-y).exp()/factorial(2*n)
    next_term = last*y/(2*n+1)
    alternating_sum = arb(1)
    term = arb(1)
    for j in range(2,42):
        term *= -y/(2*n+j)
        alternating_sum += term
    next_sum = alternating_sum - term*y/(2*n+42)
    normalized_numerator = (-2*y).exp() + next_term*alternating_sum.union(next_sum)
    upper = (1+(-2*y).exp())/2
    lower = upper - last*r**2/(1-r**2)
    assert lower > 0
    denominator = lower.union(upper)
    return normalized_numerator/denominator, last/denominator


def check_flat(n, low, high):
    lam = 2/arb(3).sqrt()
    left, right = arb(low), arb(high)
    root_sign = lambda r: mean_defect(n,r)[0] - mean_defect(n,lam*r)[0]
    assert root_sign(left) > 0 and root_sign(right) < 0
    for _ in range(35):
        mid = (left+right)/2
        value = root_sign(mid)
        if value > 0:
            left = mid
        else:
            assert value < 0
            right = mid
    root = left.union(right)
    slope_errors = []
    for area in (root, lam*root):
        defect, end_mass = mean_defect(n,area)
        # dμ/dr = 2N Var(p)/r，中间量按小缺陷重排。
        slope_error = -defect + 2*n*area*(2*defect-defect**2-end_mass)
        slope_errors.append(slope_error)
        if n == 5432:
            assert abs(slope_error) < arb("1e-2300")
            assert defect < arb("1e-2300")
    if n == 5432:
        difference = slope_errors[1]-slope_errors[0]
        assert -arb("6e-2418") < difference < -arb("5.9e-2418")
        print("四价与六价面积响应斜率之差",difference.str(8))
    print("槽位、使两面平均负载与面积成比例的六价面积区间",n,root.str(13))


def finite_distribution():
    n, x = 9, Fraction(1,5)
    weights = [(4*n*n*x)**k/factorial(2*k) for k in range(n+1)]
    for k in range(n):
        assert weights[k]*x/2 == weights[k+1]*Fraction((k+1)*(2*k+1),4*n*n)
    p = [w/sum(weights) for w in weights]
    mean = sum(Fraction(k,n)*w for k,w in enumerate(p))
    second = sum(Fraction(k*k,n*n)*w for k,w in enumerate(p))
    assert second == x*(1-p[-1])+mean/(2*n)


def verify():
    ctx.prec = 224
    finite_distribution()
    n = 5432
    theta_half = arb(1)/4
    assert theta_half.cos()>arb(".96")
    assert arb(".6")*arb(".96")-1-arb(".6").log() > arb(".08")
    error = 6*(-arb(".16")*n).exp()
    assert 2*error/(1-error) < arb("1e-370")
    # 对任一起点给 T 步的相对误差上界。偶阶链的最小邻距收缩为 1/(4N²)。
    mixing = 2*n*(-arb(n)/4 + 2*n*arb(".6")*(1-arb(".96"))).exp()
    assert mixing < arb("1e-450")
    for size in (64, 512, 5432):
        check_flat(size,".25",".28")

    cross = run_path(str(Path(__file__).with_name("ch39-regge-transverse-polarization.py")))
    base = cross["base"]
    parts = base["canonical"]()
    T,V = sp.symbols("T V")
    for cosine, expected in [(0,(T**2+T*V-8)/20),
                             (sp.Rational(1,2),(T**2-4)/20)]:
        base["canonical"] = lambda: [row for row in parts if row[1]==cosine]
        assert sp.expand(cross["cross_symbol"]()-expected) == 0


if __name__ == "__main__":
    verify()
