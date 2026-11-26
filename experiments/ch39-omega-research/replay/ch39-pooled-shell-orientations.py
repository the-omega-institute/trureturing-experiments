"""Pooled color harmonics with five uniform directions and a detectable conditional discrepancy."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

fifth = run_path(str(Path(__file__).with_name("ch39-fifth-shell-harmonic.py")))
two,q = fifth["two"],fifth["q"]


def verify():
    ctx.prec = 420
    counts60 = (1,1,19,8,31)
    counts95 = (1,1,34,10,49)
    rows = [(1,3,3),(2,4,3),(3,4,5),(4,4,5),(5,5,5)]
    count = Counter()
    for row in rows:
        for i,n in zip(row,(1,3,15)):
            count[i] += n
    assert tuple(count[i] for i in range(1,6)) == counts95
    assert min(n for n in range(56,61) if n%5==0) == 60
    scale = 1000*q(F(429,16)).exp()
    power = (two["geo"]["BOOK_PLANCK"]/scale)**2/(8*arb.pi())
    constraints = [(q(".527"),False),(q(".827"),False),(q(".827"),True),
                   (q(".9"),False),(q(".9"),True)]
    poly = fifth["constrained_polynomial"](constraints,fifth["POWERS"],-1)
    for index in range(11485,11491):
        target = two["harmonics"](acb.zeta_zero(index).imag,power,190,5)
        for counts,initial in [
            (counts60,(518,766,825,846,901)), (counts95,(485,626,824,850,901))
        ]:
            fifth["shell"]["fixed_weights"](target,counts,[q(F(v,1000)) for v in initial],10)
        moments = fifth["moment_boxes"](target,fifth["POWERS"])
        separation = sum(c*m for c,m in zip(poly,moments))/sum(abs(c) for c in poly[1:])
        print(index,"conditional harmonic discrepancy exceeds",separation.str(8))
        assert separation > q("1.33e-7")


if __name__ == "__main__":
    verify()
