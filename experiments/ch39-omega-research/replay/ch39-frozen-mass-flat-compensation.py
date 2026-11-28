"""去掉软辐边亏额时，恰保两质量的补偿与局部三质量最佳一致误差。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


def verify():
    replay = Path(__file__).resolve().parent
    model = run_path(str(replay/"ch39-two-zero-quark-window.py"))
    source = run_path(str(replay/"ch39-frozen-three-mass-local-inverse.py"))
    zero = arb(0)
    identity = arb_mat([[1,0],[0,1]])
    unit_rows = [[int(i==j) for j in range(3)] for i in range(3)]
    cases = [((i,j),arb_mat([unit_rows[i],unit_rows[j]]))
             for i,j in combinations(range(3),2)]
    cases.append(("minimax",arb_mat([[1,1,0],[0,1,1]])))
    for q in (-4,2):
        actual, evaluate = source["mass_map"](model,q)
        # 由雅可比沿线段积分包围极小质量差，避免分别包围再相减。
        def difference_at(a,b):
            _, deriv = evaluate(zero.union(a),zero.union(b),zero.union(actual))
            return deriv*arb_mat([[a],[b],[-actual]])
        full_angles = arb(0,arb("1e-12").upper())
        _, whole = evaluate(full_angles,full_angles,zero)
        for i,j in combinations(range(3),2):
            assert whole[i,0]*whole[j,1]-whole[j,0]*whole[i,1] > 0

        _, flat_j = evaluate(zero,zero,zero)
        flat_difference = difference_at(zero,zero)
        def angular(jac):
            return arb_mat([[jac[i,k] for k in range(2)] for i in range(3)])

        for label, sums in cases:
            inv = (sums*angular(flat_j)).inv()
            first = -inv*sums*flat_difference
            center = [first[i,0].mid() for i in range(2)]
            # 起点作一次 Newton 步；在邻近的有理二进小盒验证压缩自映。
            width = arb("1e-25")
            at_point = difference_at(*center)
            correction = inv*sums*at_point
            box = [v+arb(0,width.upper()) for v in center]
            assert all(abs(v) < arb("1e-12") for v in box)
            _, box_j = evaluate(*box,zero)
            derivative = identity-inv*sums*angular(box_j)
            norm = max(sum(abs(derivative[i,k]).upper() for k in range(2)) for i in range(2))
            assert max(abs(v).upper() for v in correction.entries()) + width*norm < width

            # 对包含补偿解的角盒使用均值公式，避免宽区间的质量值直接相减。
            residual = at_point+angular(box_j)*arb_mat(
                [[arb(0,width.upper())],[arb(0,width.upper())]])
            if label == "minimax":
                assert residual[0,0] < 0 and residual[1,0] > 0 and residual[2,0] < 0
                best = -residual[0,0]
                assert best > arb("1.99e-17" if q == -4 else "9.28e-17")
                print("荷、唯一局部最佳一致误差及角补偿",q,best.str(10),
                      [v.str(10) for v in box])
            else:
                other = next(i for i in range(3) if i not in label)
                assert abs(residual[other,0]) > 0
                print("荷、保持的两质量、另一质量的有符号归一偏差",q,label,
                      residual[other,0].str(9))


if __name__ == "__main__":
    verify()
