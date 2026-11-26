"""Smallest fundamental Dirac extension reversing the initial radial adjoint beta."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).parent
MATCH = run_path(str(ROOT/"ch39-quartic-soft-adjoint-match.py"))
a = MATCH["a"]


def verify():
    ctx.prec = 150
    casimir = F(12,5)
    c = 3*casimir-5
    assert c == F(11,5)
    lam, e = F(1,12), F(1,3072)
    radial = F(47,130)*e
    def yukawa_scale(n):
        return 6*casimir-F(10,3)+F(2*n,3)
    bound = yukawa_scale(3)**2/(4*c)
    gauge, box = F(375,26), F(47,65)
    assert gauge-box*bound-60*radial/F("1.16")>0
    n, active = 4, 2
    xi = yukawa_scale(n)/(c+active)
    assert xi == F(206,63) and -F(10,3)+F(2*n,3)<0
    color_moment = F(7,30)  # Tr Phi^4/(Tr Phi^2)^2 on diag(2,2,2,-3,-3)
    all_fundamentals = yukawa_scale(4)**2/(4*c)
    assert 15*color_moment+9-2*color_moment*all_fundamentals-60*color_moment*e/F("1.16")>2
    for label,kind,r,v,state in MATCH["matched_points"]():
        z = 1/state[0]
        assert z>a("1.16")
        y2 = a(xi)*z
        source = ((a(gauge)-a(box*active)*a(xi)**2)*z*z
                  +(-60+4*active*a(xi))*z*a(radial)
                  +a((F(168,25)+F(47,130)*F(64,5))*e**2)
                  +a((1-F(47,130))*30*lam**2))
        assert source < -1 and y2 < 4*arb.pi()
        # M/v has real part y*r/10. Both broken-phase masses have modulus v.
        y = y2.sqrt()
        real, imaginary = y*r/10,(1-y2*r*r/4).sqrt()
        for weight in (a("2/5"),a("-3/5")):
            assert ((real+y*r*weight)**2+imaginary**2-1).contains(0)
        print(label,kind,"16 pi^2 times radial beta",source.str(9))
    print("Fundamental Dirac copies needed in the specified one-loop class:",n)


if __name__ == "__main__":
    verify()
