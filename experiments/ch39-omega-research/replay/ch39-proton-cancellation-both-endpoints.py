"""第 177 枚正零点附近，允许八站两端同时调权时的相消边界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 192
    height = acb.zeta_zero(177).imag
    sites = (0, 2, 3, 5, 6, 9, 11, 13)
    other = tuple(i for i in sites if i != 6)
    floor = arb(85)/100
    # 各行的两个整数除以一百，给出实线性泛函的系数。
    boundaries = (
        (("-.045999", "-.045995"), 0, -1,
         ((-21,49), (-30,71), (-8,194), (-61,827), (-13,467), (-16,74), (-57,139))),
        ((".106471", ".106475"), 2, 1,
         ((-31,-73), (-20,-47), (-69,-871), (-8,-168), (-58,-141), (-16,-75), (-15,-478))),
    )
    for ends, edge, orientation, multipliers in boundaries:
        offset = arb(ends[0]).union(arb(ends[1]))
        green = shell["laplacian"](5, height + offset).inv()
        for i in sites:
            for j in sites:
                exceptional = {i, j} in ({6, 3}, {6, 5})
                real_floor = arb("-.00085") if exceptional else arb(".006")
                assert green[i, j].real > real_floor

        # 实部相消要求特殊方向的总质量 C 超过 0.87。
        assert arb(".006")/(arb(".006") + arb(".00085")) > arb(".87")
        assert arb(".13")**2/(2*arb(".87")) < arb(".01")
        reference = green[6, edge]
        assert abs(reference) > 0
        normal = orientation * reference.conjugate()/abs(reference)
        margin = arb(".00015")
        for i, (real_coefficient, imag_coefficient) in zip(other, multipliers, strict=True):
            for j in sites:
                projected = (normal * green[i, j]).imag
                for source in sites:
                    vertex = floor*green[6, j] + (1 - floor)*green[source, j]
                    bound = (real_coefficient*vertex.real + imag_coefficient*vertex.imag)/100
                    assert projected - bound > margin
        print("允许两端调权时的端点覆盖、严格投影余量", ends, margin.str(6))


if __name__ == "__main__":
    verify()
