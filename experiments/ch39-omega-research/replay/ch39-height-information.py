"""Information bounds for six homogeneous face strengths in independent finite boxes."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def covariance_bounds():
    radius, strength = F(1, 2**49), 2**115
    R2 = strength*radius**2
    sixth = 10395*108**6
    bad = sixth/R2**6
    weight = 4*1440*radius
    off = F(1, 2**30)
    assert R2/20 > 23*14+190+14*115+60
    mean = F(1092, 5)
    t = F(1, 15)
    variance = 16*36 + 288*t + 64*t*t*F(9,4)
    assert variance == F(14896,25)
    def trunc(mult, degree, moment):
        return mult*(sixth/R2**(6-degree)+moment*bad)/(1-bad)
    pointwise = 2500*strength*radius**3
    first = off+pointwise+1621*weight+trunc(15,1,108)
    second = (off+2*3242*pointwise+pointwise**2+8000000*weight
              +trunc(225,2,3*108**2))
    u = F(32, 11134**2)
    ratio = 1/(1-u)**2-1
    mixed = second + 2*mean*first + first**2
    raw_product = mean**2+9*variance+second
    two_ratios = ((1+ratio)**2-1)*(raw_product+(mean+first)**2)
    assert mixed+two_ratios < 1
    distant = F(84*14*256,2**28) + F(2*84*392,2**56)
    assert distant < F(1,100)
    central = F(192, 5)*variance
    assert 22840 < central-34-distant <= central+34+distant < 23000


def classify():
    ctx.prec = 160
    setup = run_path(str(Path(__file__).with_name("ch39-full-color-obstruction.py")))
    q = setup["q"]
    power = (setup["high"]["geo"]["BOOK_PLANCK"]/(1000*q(setup["L0"]).exp()))**2/(8*arb.pi())
    strengths = [1/(acb.zeta_zero(i).imag**2+q(1)/4) for i in range(11485,11491)]
    gaps = [power*(a-b) for a,b in zip(strengths,strengths[1:])]
    assert min(gaps) > q("0.0000204") and gaps[-1] < q("0.00002042")
    expected_gap = 22840*F(204,10**7)
    assert expected_gap > F(465,1000)
    # 1/1000 of kappa^-1 calibration allowance on either midpoint.
    assert F(460,180000*(F(465,2000)-F(1,1000))**2) < F(1,20)
    assert (F(9,5)**2)/(50*23000*F(2042,10**8)**2) > 6700
    print("95% discrimination: necessary sites times boxes > 6700*kappa^2; sufficient 180000*kappa^2.")


if __name__ == "__main__":
    covariance_bounds()
    classify()
