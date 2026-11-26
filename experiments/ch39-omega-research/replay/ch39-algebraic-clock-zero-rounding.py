"""从第 36、1 枚实际高度用同一个最大分母确定有理钟面。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    graph = run_path(str(root / "ch39-arc-relative-determinant.py"))
    ctx.prec = 256
    choices = sorted({Fraction(a, b) for b in range(1, 30) for a in range(-b, b+1)})
    targets = iter(((-1, 0, 1), (-1, 1, 29), (1, -4, 17), (1, 1, 6)))
    margins = []
    for n in (36, 1):
        phases = graph["phases"](acb.zeta_zero(n).imag)
        for d in (724, 804):
            sign, a, b = next(targets)
            z = phases[d]*sign
            assert z.real > 0
            h = z.imag/(1+z.real)
            nearest = Fraction(a, b)
            dist = abs(h-arb(a)/b)
            other = min(abs(h-arb(v.numerator)/v.denominator).lower() for v in choices if v != nearest)
            margin = (other-dist)/2
            assert margin > 0
            margins.append(margin)
            if nearest == Fraction(-4, 17):
                assert abs(h+arb(7)/30) < dist
            print("编号、弦长、半角最邻近有理数、到任何舍入分界的距离下界",
                  n, d, str(nearest), margin.lower().str(12))
    assert min(v.lower() for v in margins) > arb("1e-6")


if __name__ == "__main__":
    verify()
