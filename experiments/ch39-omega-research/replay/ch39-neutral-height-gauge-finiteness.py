"""用整数行列式和 Sturm 计数核验中性质量首次退化的位置。"""
from pathlib import Path
from runpy import run_path

from flint import fmpz_mat
import sympy as sp


def verify():
    shell = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))
    fixed, slope = fmpz_mat(19, 19), fmpz_mat(19, 19)
    for i, j, _, _ in shell["geo"]["edges"]:
        spoke = shell["soft"]["pole"] in (i, j) and (i in shell["soft"]["star"]["orbit"]
                                                    or j in shell["soft"]["star"]["orbit"])
        for matrix, weight in ((fixed, 3), (slope, -3*int(spoke))):
            matrix[i, i] += weight
            matrix[j, j] += weight
            matrix[i, j] -= weight
            matrix[j, i] -= weight
    for i in range(19):
        fixed[i, i] -= 45
        slope[i, i] += 1
    t = sp.symbols("t")
    quartic = 2*t**4-3*t**3-333*t**2-1998*t-1296
    other = 11*t**4-156*t**3+234*t**2+2376*t-5184
    determinant = -(t-45)*(t+6)**3*(t+12)**3*quartic**2*other
    # 两边次数至多 19，二十个互异整数点即给出多项式恒等式。
    for k in range(20):
        assert int((fixed+k*slope).det()) == determinant.subs(t, k)
    left, right = sp.Rational("-0.738688949"), sp.Rational("-0.738688948")
    assert sp.Poly(quartic, t).count_roots(left, right) == 1
    assert sp.Poly(quartic, t).count_roots(right, 0) == 0
    assert sp.Poly(other, t).count_roots(left, 0) == 0
    assert determinant.subs(t, 0) != 0
    beta = (sp.Rational(41, 10), -sp.Rational(19, 6), -7)
    r = (beta[0]-beta[2])/(beta[0]-beta[1])
    projection = -7*r+sp.Rational(21, 2)
    assert projection == -sp.Rational(21, 109)
    print("最近负根所在区间", left, right, "中性重矢量的规范投影", projection)


if __name__ == "__main__":
    verify()
