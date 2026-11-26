"""Indices and uniform local beta-function signs for the full matched product group."""
from fractions import Fraction as F
from flint import arb, ctx, fmpq


def coefficients(ca, fermions, scalars):
    b = -F(11,3)*ca+F(2,3)*sum(t for t,_ in fermions)+F(1,3)*sum(t for t,_ in scalars)
    diagonal = (-F(34,3)*ca**2+sum(t*(2*c+F(10,3)*ca) for t,c in fermions)
                +sum(t*(4*c+F(2,3)*ca) for t,c in scalars))
    return b, diagonal


def target_inverse_coefficients():
    # X and Y each supply ten fundamental complex scalars to SU(10).
    index = 2*10*F(1,2)
    one, diagonal = coefficients(10, [(3,F(99,20))], [(index,F(99,20))])
    cross = 4*index*F(18,5)
    return tuple(fmpq(c.numerator,c.denominator) for c in (-one,-diagonal,cross))


def verify():
    ctx.prec = 140
    five = coefficients(5,[(F(3,2),F(12,5)),(F(9,2),F(18,5))],
                        [(30,F(18,5)),(F(1,2),F(12,5))])
    ten, diagonal, cross = target_inverse_coefficients()
    assert five == (F(-25,6),F(5921,15)) and (ten,diagonal,cross) == (fmpq(94,3),fmpq(22169,30),144)
    assert 4*30*F(99,20) == 594 and 24*594 == 99*cross
    assert 100*F(18,5)/24 == 15 and 100*F(99,20)/99 == 5
    pi = arb.pi()
    cap = (4*pi).sqrt()
    e = arb("4.61e-5")
    d = 1000/arb("7.63e14")
    k = (e*cap)**3/d
    assert k < cap**2
    ysum = 4*cap**2+(k/cap)**2+(e*cap)**2
    gauge = arb(".618")**2*((arb(5921)/15).sqrt()+arb(3*594).sqrt())**2
    reduced = -arb(25)/6+(gauge-15*ysum)/(16*pi**2)
    assert ysum < arb("51.15") and reduced > arb(".296")
    assert -ten+cross*4*pi/(16*pi**2) < -19
    print("b=(-25/6,-94/3); B=((5921/15,594),(144,-22169/30)).")
    print("Uniform lower bound on the reduced SU(5) beta",reduced.str(11))


if __name__ == "__main__":
    verify()
