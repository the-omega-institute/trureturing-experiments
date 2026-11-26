"""Uniform lower and reweighting bounds for the finite three-family insertion."""
from fractions import Fraction as F
from itertools import product
from flint import acb, arb, ctx


def verify():
    ctx.prec = 140
    t = F(1, 2**20)
    change = 32*t
    assert F(1,100)-change > F(9,1000) and change/F(9,1000) < F(1,250)
    # Last certificate: massive D^dagger D eigenvalues > 0.018 - 0.0041.
    assert F(18,1000)-F(41,10000) > (F(1,10)+F(1,250))**2
    assert 4*16-(16-1) == 49
    faces = 16*sum(1 for x in product((0,1,2),repeat=4) if 1 in x and 2 in x)
    assert faces == 800
    c = arb(30)+arb(2)/15
    strengths = []
    for j in range(11485,11491):
        z = acb.zeta_zero(j)
        assert z.real == arb(1)/2
        strengths.append(1/(z.imag**2+arb(1)/4))
    smax, smin = strengths[0], strengths[-1]
    r = (1-smax*c)/(1+smax*c)
    assert r**faces > arb(1)/2
    bound = F(216,2*16**5)*(t/4)**(3*49)*F(1,20)**1440
    assert bound > F(1,2**9472)
    delta = 2*faces*c*(smax-smin)/(1-c*smax)
    turn_on = 2*faces*c*smax/(1-c*smax)
    assert delta.exp()-1 < arb("3e-7") and turn_on.exp()-1 < arb(".0004")
    print("2^4 normalized RH-weighted insertion > 2^-9472.")
    print("Max six-label relative reweighting bound", (delta.exp()-1).str(12))


if __name__ == "__main__":
    verify()
