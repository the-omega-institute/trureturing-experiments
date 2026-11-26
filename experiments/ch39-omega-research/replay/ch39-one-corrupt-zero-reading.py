"""共同零位未知且任意一个零点的读数可能出错时的锐前缀。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def verify():
    source = Path(__file__).with_name("ch39-winding-difference-prefix.py")
    prefix = run_path(str(source))
    ctx.prec = 320
    equalities, action, _, _ = prefix["rows"]()
    fit = arb_mat([row + [0] for row in equalities] +
                  [row + [1] for row in action[:101]])
    prediction = arb_mat([row + [1] for row in action[101:103]]) * fit.inv()

    # 末两项各自减去由前 101 项与五条精确等式给出的预测。
    # 每个观测位置的单一错误，都沿以下矩阵的一列改变这两个残差。
    columns = [(-prediction[0, 5+j], -prediction[1, 5+j]) for j in range(101)]
    columns.extend([(arb(1), arb(0)), (arb(0), arb(1))])
    norms = [(a*a + b*b).sqrt() for a, b in columns]
    assert all(norm > 0 for norm in norms)
    sines = []
    for j, k in combinations(range(len(columns)), 2):
        a, b = columns[j]
        c, d = columns[k]
        sine = abs(a*d - b*c)/(norms[j]*norms[k])
        assert sine > 0, (j+1, k+1)
        sines.append((sine, j+1, k+1))
    smallest, j, k = min(sines, key=lambda entry: float(entry[0]))
    assert all(j == u and k == v or smallest < sine for sine, u, v in sines)
    assert len(sines) == 5253
    print("最接近的误差方向对应编号", j, k, "锐角正弦", smallest.str(20))
    print("校验方向列范数范围", min(norms).str(14), max(norms).str(14))


if __name__ == "__main__":
    verify()
