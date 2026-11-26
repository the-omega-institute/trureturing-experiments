"""Flat discrete holonomy and the mixed finite Weyl index."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path
from flint import arb, arb_poly, ctx

flux = run_path(str(Path(__file__).with_name("ch39-local-triality-mixed-flux.py")))
ROWS = {
    -1: (0,1,0,1,0, 1,1,1,0,1, 1,0,1,1,1, 1,1,1,1,2, 1,0,0,1,0),
     1: (0,1,1,1,1, 1,0,1,1,0, 1,1,0,1,1, 1,1,1,0,1, 1,0,1,1,1),
}


def verify():
    ctx.prec = 640
    cut = arb(1)/250
    values = {}
    for twist, expected in ROWS.items():
        pattern = tuple((a+twist) % 3 for a in flux["FIRST"])
        indices = []
        for matrix in flux["phase_fibers"](pattern, flux["SECOND"]):
            poly = arb_poly([c.real for c in matrix.charpoly().coeffs()])
            count = flux["base"]["positive_roots"](poly, cut)
            assert count == flux["base"]["positive_roots"](poly, -cut)
            indices.append(50-count)
        assert tuple(indices) == expected
        values[twist] = sum(indices)
    assert values == {-1:18, 1:19} and 2*(values[1]-values[-1]) == 2
    assert Fraction(256, 3*2**14) < Fraction(1,100)
    print("Both flat twists certified at gap > 1/250; indices 18 and 19, triality exponent 2.")


if __name__ == "__main__":
    verify()
