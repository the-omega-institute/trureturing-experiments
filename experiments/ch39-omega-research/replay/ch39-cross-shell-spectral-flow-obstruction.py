"""计算荷六零弱差所需的端点谱相位之外的酉变化与一次交换的时间积下界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adjoint(a):
    return a.conjugate().transpose()


def verify(root):
    records, _ = run_path(str(root / "ch39-three-shell-two-space-readout.py"))["verify"](root, return_data=True)
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    spec = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    center = 2_306_888_100_085
    ctx.prec = 320
    cutoff = 16*arb("91.1876")/arb("1.22089e19")
    for side in (-1, 1):
        a, b = center + side, center + 2*side
        energy, source = spec["spectrum"](shell, 6, a)
        _, target = spec["spectrum"](shell, 6, b)
        weights = spec["training_values"](energy, records)
        high = [j for j, w in enumerate(weights) if w > cutoff]
        assert len(high) == (4 if side == -1 else 2)
        m = min(w.lower() for w in weights)
        h = min(weights[j].lower() for j in high)
        cap = ((cutoff-m)/(h-m)).sqrt()
        line = acb_mat([[target[j, 15]] for j in range(19)])
        overlap = adjoint(source)*line
        before = sum(abs(overlap[j, 0])**2 for j in high).sqrt()
        # Δ=角(目标第16线,源低空间)−允许的最大角。
        angle = before.asin()-cap.asin()
        distance = 2*(angle/2).sin()
        A, B = [shell["laplacian"](6, acb.zeta_zero(n).imag) for n in (a,b)]
        commutator = A*B-B*A
        # 对半正定 C* C，最大本征值的四次方不超过四次幂的迹。
        square = adjoint(commutator)*commutator
        fourth = square*square
        size = (fourth*fourth).trace().real.root(8).upper()
        assert distance > 0 and size > 0
        print("编号、必需酉距离下界、主角差、交换子范数上界、两端时间积下界",
              a, b, distance.str(13), angle.str(13), size.str(13),
              (distance/size).str(13), flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
