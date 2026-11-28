"""首枚正零点上，一端沿极点辐边同相制备时，另一端精确相消的最小相位宽度。"""
from itertools import combinations_with_replacement, product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 192
    green = shell["laplacian"](5, acb.zeta_zero(1).imag).inv()
    sites = (0, 2, 3, 5, 6, 9, 11, 13)
    for i, j in product(sites, repeat=2):
        assert green[i, j].real > 0 and abs(green[i, j].arg()) < arb(".82")

    ref = green[6, 5].conjugate()*green[6, 0]
    width = ref.arg()
    assert arb(".84") < width < arb(".85")
    for left, right in product(sites, repeat=2):
        def coefficient(i, m):
            mixed = green[i, left].conjugate()*green[m, right]
            mixed += green[m, left].conjugate()*green[i, right]
            return -(ref.conjugate()*mixed).imag/2

        for i, m in combinations_with_replacement(sites, 2):
            if (left, right, i, m) == (5, 0, 6, 6):
                # 此项为 -Im(conj(ref)*ref)，恒为零。
                continue
            if (left, right, i, m) == (9, 6, 3, 5):
                diag = coefficient(3, 3)*coefficient(5, 5)
                assert coefficient(3, 5)**2 < diag
                continue
            assert coefficient(i, m) > 0

    threshold = arb.pi()-width
    assert arb("2.2984473859") < threshold < arb("2.2984473861")
    print("接收相位最小宽度（弧度）", threshold.str(15))


if __name__ == "__main__":
    verify()
