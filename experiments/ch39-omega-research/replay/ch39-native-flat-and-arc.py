"""中心奇偶一圈模型的平坦对照与两种有限整数读数检验。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

model = run_path(str(Path(__file__).with_name("ch39-native-center-parity-crossing.py")))


def projection(values):
    y, w, c = values
    return (115*y - 333*w + 218*c)/218


def verify():
    ctx.prec = 256
    original, radius = model["setup"]()
    flat = {charge: original[0] for charge in original}
    mz, beta = model["mz"], model["beta"]
    min_log = (arb("1e13")/mz).log()

    def intercept(levels, u):
        return model["intercept"](levels, radius, u, min_log)

    # Endpoint signs hold uniformly over the full finite subtraction window.
    def bracket(levels, lower, upper):
        lo, hi = arb(lower), arb(upper)
        assert projection(intercept(levels, lo)) < 0 < projection(intercept(levels, hi))
        return lo, hi

    lo, hi = bracket(flat, "2.85559750565", "2.85559750569")
    y, w, _ = intercept(flat, lo.union(hi))
    scale_log = 2*arb.pi()*(y - w)/(beta[0] - beta[1])
    scale = mz*scale_log.exp()
    common = y - beta[0]*scale_log/(2*arb.pi())
    assert arb("6.04924187e13") < scale < arb("6.04924189e13")
    assert arb("38.18972539") < common < arb("38.18972541")
    print("平坦相位的刚度、交汇尺度与共同逆耦合", lo.union(hi), scale, common)

    lo, hi = bracket(original, "2.83087012440955", "2.83087012441002")
    fitted = lo.union(hi)
    residual = projection(intercept(flat, fitted))
    assert arb("-.046985") < residual < arb("-.046983")
    print("原刚度置于平坦连接的投影", residual)

    # A real eigenvalue of absolute value u>1 in GL(2,Z) has
    # |trace| = u+1/u or u-1/u, according to the determinant.
    for trace in (fitted + 1/fitted, fitted - 1/fitted):
        assert not trace.contains_integer()
        print("整系数幺模矩阵所需的迹绝对值", trace)

    # cot(phi/2) fixes cos(phi) monotonically for the fitted positive u.
    def cosine(u):
        return (u*u - 1)/(u*u + 1)

    lower = Fraction("0.77811862852057")
    upper = Fraction("0.77811862852065")
    assert arb(lower.numerator)/lower.denominator < cosine(lo)
    assert cosine(hi) < arb(upper.numerator)/upper.denominator
    chord = 1734*(1 - cosine(lo).union(cosine(hi)))
    assert 384 < chord < 385
    # With pole (17D,17D,17D), S=X+Y+Z=51D cos(phi).
    # S and D have the same parity because S = sum X_i^2 mod 2.
    # (S-D)/(2D) belongs to this Farey interval of determinant one.
    left, right = Fraction(2252050, 116433), Fraction(19659447, 1016411)
    assert left < (51*lower-1)/2 < (51*upper-1)/2 < right
    assert right.numerator*left.denominator-left.numerator*right.denominator == 1
    min_denominator = left.denominator + right.denominator
    assert min_denominator == 1132844
    # If D is even, the sum of three squares is 0 mod 4, so every
    # coordinate is even. Halving the coordinates and D preserves the angle.
    assert all(sum(row) % 4 for row in ((1,0,0),(1,1,0),(1,1,1)))
    print("原壳平方弦长、可能的最小整数放大倍数下界", chord, min_denominator+1)


if __name__ == "__main__":
    verify()
