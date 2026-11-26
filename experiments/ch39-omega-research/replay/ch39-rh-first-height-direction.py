"""Three offset readings force three real probe directions at the first true RH height."""
from pathlib import Path
from runpy import run_path

from flint import acb,arb,ctx

directory=Path(__file__).parent
read=run_path(str(directory/"ch39-finite-nbb-coherent-four-windows.py"))["readings"]


def main():
    ctx.prec,ctx.cap=256,3
    f=read(3)
    assert acb.zeta_zero(1).imag<15 and acb.zeta_zero(2).imag>21
    g=(3*f[0]-4*f[1]+f[2])/12
    z=arb(1)/21**2
    ceiling=z*z/((1+z)*(1+4*z))
    margin=(g-ceiling)*12/5
    assert margin>arb("0.00000243644")
    print("At least one operator error exceeds 2.43644e-6 if the first height has rank below three.")
    # With zero momentum normalized, a direction whose first-height weight is
    # tiny must cover the same shortfall. The first atom sets the largest filter.
    zfirst=1/acb.zeta_zero(1).imag**2
    znext=1/acb.zeta_zero(2).imag**2
    val=lambda x:x*x/((1+x)*(1+4*x))
    weight=(g-val(znext))/(val(zfirst)-val(znext))
    print("Minimum first-height Gram eigenvalue from these three exact readings",weight.str(14))
    # At the first two offsets alone, one rank-one first-height residue suffices:
    # the other two axes share the second low height, and all three axes
    # use separately chosen high heights.
    target=f[1]
    lows=[acb.zeta_zero(j).imag**2 for j in (1,2)]
    highs=[acb.zeta_zero(j).imag**2 for j in (5,6,7)]
    low_response=[g/(g+1) for g in lows]
    assert max(low_response)<target<min(g/(g+1) for g in highs)
    print("Offsets 0,1 alone allow three axes using ranks one and two at the first two heights.")


if __name__=="__main__":
    main()
