"""四维共享边的 Haar 刚度、稀释极限与有限面强度界。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    ctx.prec = 160
    reps = run_path(str(Path(__file__).with_name("ch39-native-trace-plaquette-jets.py")))
    family, index, dim = (reps[k] for k in ("family", "index", "dimension"))
    bare = [sum(dim(r) * index(r)[i] for r in family) for i in range(3)]
    schur = [sum(index(r)[i] for r in family) / Fraction(8) for i in range(3)]
    assert bare == [120, 2, 2] and schur == [Fraction(33, 4), Fraction(1, 16), Fraction(1, 16)]
    # All representations and their conjugates are inequivalent, already by hypercharge.
    assert len({sign*r[0] for r in family for sign in (1, -1)}) == 2*len(family)
    assert sum(dim(r) for r in family) == 15
    d_y = sum(dim(r) * abs(r[0]) for r in family) / 2
    d_w = sum(dim(r) for r in family if r[1] == "two") / 4
    assert (d_y, d_w) == (18, 2)
    # Divide L² norms by tau = tr(F_q), which is common across the five reps.
    y, w = (arb(a.numerator).sqrt() / arb(a.denominator).sqrt() for a in schur[:2])
    for tmax, emax, bound in [(19, 2500, "65.84"), ("7.311", 1000, "66.41"),
                              (3, 1000, "97.91")]:
        epsilon = arb(1) / emax
        error = epsilon / (1-epsilon)
        ratio = (1-epsilon)**(90*arb(tmax)) * ((y-d_y*error)/(w+d_w*error))**2
        assert y > d_y*error and w > d_w*error and ratio > arb(bound)
        print("共同迹上界、面强度上限、共享边超荷/弱刚度比的下界",
              tmax, "1/" + str(emax), ratio.str(12))


if __name__ == "__main__":
    verify()
