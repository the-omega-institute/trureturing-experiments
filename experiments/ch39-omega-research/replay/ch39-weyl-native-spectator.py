"""用同一交点能量计数第二份十九站壳的局部第二陈数。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    cone = run_path(str(root / "ch39-first-chern-edge.py"))
    cone["verify"](root)
    native = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    levels = run_path(str(root / "ch39-weyl-light-yukawa-determinant.py"))["spectator_levels"]
    ctx.prec = 256
    center = arb("80.928425691033686017")
    weights = [arb(".01615054064436668415"), arb("11.16447955978410845")]
    width = arb("1e-12")
    energy = arb("8.929314082054833696") + arb(0, width)
    clock_rates = [abs(native["geo"]["theta"](d)-native["geo"]["theta"](244))
                   for _, _, d, _ in native["geo"]["edges"] if d in (724, 804)]
    for q, count in ((6, 3), (-4, 3), (2, 4), (1, 4)):
        eig = levels(native, cone["EDGES"], center, weights, q)
        shift = width*(abs(q)*sum(clock_rates) + 5)
        below = eig[count-1].upper() + shift
        above = eig[count].lower() - shift
        assert eig[0] > shift and below < energy < above
        if q == 6:
            assert below/energy < arb(".983") and above/energy > arb("1.0125")
        print("第二份壳的荷、球面第二陈数、包围半径除以交点能量的可行开区间",
              q, count, (below/energy).upper().str(14), (above/energy).lower().str(14))


if __name__ == "__main__":
    verify()
