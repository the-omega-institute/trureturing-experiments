"""用两个相对谱迹确定任意共同复源的九项与十项最小最大误差。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def normalized_traces(up, down):
    relative = up * down.inv()
    return tuple(sum(abs(z)**2 for z in matrix.entries()) / abs(matrix.det())**(arb(2)/3)
                 for matrix in (relative, relative.inv()))


def verify():
    replay = Path(__file__).resolve().parent
    library = run_path(str(replay / "ch39-nine-window-cp-cosine.py"))
    affine_trace = library["affine_trace"]
    dual = library["Dual"]
    model, up_box, down_box = run_path(str(replay / "ch39-two-zero-orbit-condition.py"))["calibration"]()
    at_zero = run_path(str(replay / "ch39-quadratic-gram-source.py"))["at_zero"]
    centers, widths = model["ref"][:9], model["tol"][:9]
    signs = (-1, 1, 1, -1, 1, 1, 1, 1, 1)

    def corner(radius):
        return [center + sign*radius*width
                for center, width, sign in zip(centers, widths, signs, strict=True)]

    def residual(values, direct, reciprocal):
        base, slope = affine_trace(values)
        inverse_base, inverse_slope = affine_trace(values, inverse=True)
        return inverse_base + inverse_slope*((-base+direct)/slope) - reciprocal

    traces = normalized_traces(up_box, down_box)
    box = [c + arb(0, (arb(".059")*w).upper()) for c,w in zip(centers,widths,strict=True)]
    for j, sign in enumerate(signs):
        values = [dual(v, int(i==j)) for i,v in enumerate(box)]
        assert sign*residual(values, *traces).d > 0

    actual = normalized_traces(at_zero(model, -4, 2345666), at_zero(model, 2, 964890))
    assert all(wide.contains(narrow) for wide,narrow in zip(traces,actual,strict=True))
    for label, pair, (lo,hi) in (
            ("原制备盒", traces, (".05739", ".05748")),
            ("实际零点", actual, (".0574339450674740", ".0574339450674742"))):
        left,right = arb(lo),arb(hi)
        assert residual(corner(left), *pair) < 0 < residual(corner(right), *pair)
        cosine = library["phase"](corner(left.union(right)), pair[0])
        assert -1 < cosine < 1
        delta = cosine.acos()
        assert abs(delta-model["ref"][9]) < left*model["tol"][9]
        if label == "实际零点":
            assert arb("1.13891617544348") < delta < arb("1.13891617544349")
        print(label, "最优半宽倍数", lo, hi, "正 CP 相位", delta.str(18))


if __name__ == "__main__":
    verify()
