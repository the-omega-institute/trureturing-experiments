"""按整数编号在反自伴 Cayley 坐标上外推一步。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def extrapolate(x,y):
    eye = acb_mat([[int(i == j) for j in range(19)] for i in range(19)])
    a, b = eye+x, eye+y
    assert not a.det().contains(0) and not b.det().contains(0)
    # 2(I+Y)^(-1)-(I+X)^(-1) = b^(-1)(2a-b)a^(-1).
    # 先消去两个图表的逆，减少区间依赖造成的包围膨胀。
    denominator = 2*a-b
    assert not denominator.det().contains(0)
    inverse = denominator.mid().inv()
    residual = eye-inverse*denominator
    error = sum(abs(z).upper()**2 for z in residual.entries()).sqrt()
    assert error < 1
    # M^{-1}b=(I-E)^{-1}(Cb)，量化右端而非先包围整个逆。
    solved = inverse*b
    defect = inverse*(b-denominator*solved)
    correction = sum(abs(z).upper()**2 for z in defect.entries()).sqrt()/(1-error)
    radius = arb(0,correction.upper())
    solved = acb_mat([[solved[i,j]+acb(radius,radius) for j in range(19)] for i in range(19)])
    step = a*solved-eye
    print("预条件残差、修正半径、第三步最大元素半径",
          error.str(8), correction.str(8),max(z.rad().upper() for z in step.entries()).str(8),flush=True)
    assert (step.conjugate().transpose()*step).contains(eye)
    return step


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    run_path(str(root/"ch39-fourth-arithmetic-bridge-holdout.py"))["verify"](root,extrapolate)
