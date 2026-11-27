"""第六壳两个固定重耦合截面的 CP 行列式可达与排除。"""

from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def adjoint(a):
    return a.conjugate().transpose()


def direct_compression(v, h, a, b, blocks):
    """以法向量的正交补为列空间，独立构造两枚三阶矩阵。"""
    leftover = (1 - h).sqrt()
    frame = acb_mat(
        [
            [-v[1].conjugate(), leftover * v[0], 0],
            [v[0].conjugate(), leftover * v[1], 0],
            [0, -(h * a).sqrt() / leftover, b.sqrt() / leftover],
            [0, -(h * b).sqrt() / leftover, -a.sqrt() / leftover],
        ]
    )
    matrices = []
    for trivial, first, second in blocks:
        full = acb_mat(
            [
                [trivial[0, 0], trivial[0, 1], 0, 0],
                [trivial[1, 0], trivial[1, 1], 0, 0],
                [0, 0, first, 0],
                [0, 0, 0, second],
            ]
        )
        matrices.append(adjoint(frame) * full * frame)
    return matrices


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    lib("common-plane-cp-factor")["verify"]()
    _, target = lib("sixth-cp-commutator-amplitude")["verify"](root)
    (tu, td), ((ua, ub), (da, db)) = lib("sixth-heavy-frontier")["blocks"](root)
    x = arb(".8659")
    assert ua < x < ub and db < da
    eye = acb_mat([[1, 0], [0, 1]])
    comm = tu * td - td * tu
    norm = comm.det().real.sqrt()
    det = ((tu - eye * ua) * (da - db) + (td - eye * da) * (ub - ua)).det()
    assert det.imag.contains(0) and det.real > 0 and norm > 0
    alpha = (x - ua) / (ub - ua)
    ru = (eye * x - tu).inv() * (ub - ua) * alpha * (1 - alpha)
    for y, reachable in ((arb(".01653"), True), (arb(".016525"), False)):
        assert db < y < da
        beta = (y - db) / (da - db)
        sigma = alpha + beta - 1
        assert sigma > 0
        rd = (eye * y - td).inv() * (da - db) * beta * (1 - beta)
        total = ru + rd
        # 二阶 Hermitian 矩阵最小本征值，使用对角差与非对角模的平方和。
        center = total.trace().real / 2
        radius = (
            ((total[0, 0] - total[1, 1]).real / 2) ** 2 + abs(total[0, 1]) ** 2
        ).sqrt()
        smallest = center - radius
        assert smallest > 0
        h_upper = sigma / (sigma + smallest)
        absolute_bound = (h_upper * det.real * norm / 4).upper()
        print(
            "底耦合、斜边余量、平凡压缩核最小本征值下界、CP 绝对值上界",
            y.str(8),
            sigma.lower().str(12),
            smallest.lower().str(12),
            absolute_bound.str(12),
            flush=True,
        )
        if not reachable:
            assert absolute_bound < arb("4.042e-15") < target
            continue
        for sign in (-1, 1):
            # i[T_u,T_d] 的两个本征射线给出相反的 CP 符号。
            projector = (eye + comm * (acb(0, sign) / norm)) / 2
            su = (projector * ru).trace().real
            sd = (projector * rd).trace().real
            denom = sigma + su + sd
            h = sigma / denom
            a = ((1 - beta) * su + alpha * sd) / denom
            b = (beta * su + (1 - alpha) * sd) / denom
            assert min(h, a, b) > 0 and (h + a + b).contains(1)
            v = [projector[i, 0] / projector[0, 0].real.sqrt() for i in range(2)]
            v_column = acb_mat([[v[0]], [v[1]]])
            assert (adjoint(v_column) * v_column)[0, 0].real.contains(1)
            assert (
                a / (x - ua) + b / (x - ub) + h * su / ((ub - ua) * alpha * (1 - alpha))
            ).contains(0)
            assert (
                a / (y - da) + b / (y - db) + h * sd / ((da - db) * beta * (1 - beta))
            ).contains(0)
            predicted = sign * a * b * h * det.real * norm
            assert sign * predicted > arb("2.8e-14") > target
            left, right = direct_compression(v, h, a, b, ((tu, ua, ub), (td, da, db)))
            identity3 = acb_mat([[int(i == j) for j in range(3)] for i in range(3)])
            for eigen, matrix in ((x, left), (y, right)):
                poly = (identity3 * eigen - matrix).det()
                assert poly.real.contains(0) and poly.imag.contains(0)
            independent = (left * right - right * left).det()
            assert independent.real.contains(0) and independent.imag.overlaps(predicted)
            assert sign * independent.imag > arb("2.90e-14")
            print(
                "方向、因子式读数、直接压缩区间",
                sign,
                predicted.lower().str(13),
                predicted.upper().str(13),
                independent.imag.lower().str(13),
                independent.imag.upper().str(13),
                flush=True,
            )


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
