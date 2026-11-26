"""三枚相邻实际零点共用的普朗克有界截止所需负质量的双侧证书。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat


def verify(root):
    fields = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    n = 2_306_888_100_085
    single, lib, spectra = fields(root, (n, n + 1, n + 2))
    bounds = run_path(str(root / "ch39-common-cutoff-negative-mass.py"))
    epsilon = arb("1e-122")
    target = arb_mat([[1], [1], [epsilon], [0], [0], [0], [0]])

    def column(x):
        readings = [single["moments"](fields, lib["expint"], x) for fields in spectra]
        k, w = readings[0]
        return [arb(1), k, w, readings[1][0]-k, readings[1][1]-w,
                readings[2][0]-k, readings[2][1]-w]

    sites = list(map(arb, ("7.416959887e-10", "1.45768048386e-7", "8.31773636916e-7",
                           "1.90167299418e-5", ".0332693440288", ".547496172832", "1")))
    basis = arb_mat([column(x) for x in sites]).transpose()
    assert not basis.det().contains(0)

    # 用中点矩阵的近似逆预处理，并以原方程的区间残差严格包住权重校正。
    def midpoints(mat):
        return arb_mat([[mat[i, j].mid() for j in range(mat.ncols())] for i in range(mat.nrows())])

    pre = midpoints(midpoints(basis).inv())
    assert not pre.det().contains(0)
    initial = midpoints(pre*target)
    weights = (initial + (pre*basis).solve(pre*(target-basis*initial))).entries()
    signs = (-1, 1, -1, 1, -1, 1, -1)
    assert all(sign*w > 0 for sign, w in zip(signs, weights, strict=True))
    upper = -sum(weights[::2])

    # 对偶系数可以任取。插值解取精确二进制中点后，在所有端点和转折窗重新计算函数值。
    dual = midpoints(basis.transpose().solve(arb_mat([[-1], [0], [-1], [0], [-1], [0], [-1]]))).entries()
    _, alpha, beta, gamma1, eta1, gamma2, eta2 = dual
    multipliers = ((alpha-gamma1-gamma2, beta-eta1-eta2), (gamma1, eta1), (gamma2, eta2))
    rows = [(z, a*k/(2*arb.pi()), -g*w/(4*arb.pi()))
            for fields, (k, w) in zip(spectra, multipliers, strict=True) for g, a, z in fields]

    def face(x):
        return sum(c*v for c, v in zip(dual, column(x), strict=True))

    def small(x):
        return bounds["small_derivative"](x, rows)

    large = bounds["long_series"](rows, 27)
    cut = arb("1e-5")
    centers = [arb("4.57437940828e-11"), *sites[:-1]]
    radii = list(map(arb, ("1e-13", "1e-16", "2e-15", "1e-14", "4e-13", "2e-10", "1e-9")))
    left, sign = arb(0), 1
    ranges = [face(left), face(arb(1))]
    leaves = []
    for x, radius in zip(centers, radii, strict=True):
        right = x-radius
        derivative = small if x < cut else large
        if left < cut < right:
            leaves.append(bounds["monotone"](left, cut, sign, small))
            left = cut
        leaves.append(bounds["monotone"](left, right, sign, derivative))
        left = x+radius
        error = radius*abs(derivative(right.union(left)))
        ranges.append(face(x) + arb(0, error.upper()))
        sign = -sign
    leaves.append(bounds["monotone"](left, arb(1), sign, large))
    bottom = min(v.lower() for v in ranges)
    top = max(v.upper() for v in ranges)
    lower = (dual[0]+alpha+epsilon*beta-top)/(top-bottom)

    assert arb("2911897.02372933") < lower
    assert upper < arb("2911897.02372935")
    ratio_lo = lower/arb("337760.30905793")
    ratio_hi = upper/arb("337760.30905784")
    assert arb("8.62119362648") < ratio_lo and ratio_hi < arb("8.62119362650")
    print("七个原子有符号权", *[w.str(20) for w in weights])
    print("各单调段的细分数、对偶全局值域", leaves, bottom.str(25), top.str(25))
    print("负质量的严格下界与可行上界", lower.str(23), upper.str(23))
    print("三高度与双高度的下确界比", ratio_lo.str(18), ratio_hi.str(18))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
