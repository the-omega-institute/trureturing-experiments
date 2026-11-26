"""Minimum uniform cardinalities for six zero-weight harmonics on a finite Behrend shell."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb

prior = run_path(str(Path(__file__).with_name("ch39-higgs-sixth-harmonic.py")))
q, arc, powers = prior["q"], prior["arc"], prior["powers"]
CUTS = (
    (488462,811783,892422,949366), (488412,811761,892410,949352),
    (488335,811728,892392,949330), (488302,811713,892384,949321),
    (488233,811682,892367,949301), (488200,811668,892359,949291),
)


def verify():
    two = arc["fifth"]["two"]
    interpolate = arc["fifth"]["constrained_polynomial"]
    double, value, mean = prior["double"], prior["polynomial"], prior["expected"]
    scale = 1000*q(F(429,16)).exp()
    power = (two["geo"]["BOOK_PLANCK"]/scale)**2/(8*arb.pi())
    sizes = (1,3,3,6,6)
    sectors = ((0,1,2,2,2),(5,4,4,2,4),(3,4,4,4,4))
    counts = [0]*6
    for rank,row in zip((1,1,2),sectors):
        for size,col in zip(sizes,row):
            counts[col] += rank*size
    assert counts == [1,3,21,2,48,1]
    for zero,cuts in enumerate(CUTS,11485):
        a,b,c,h = cuts
        A,B,C,D,H = [q(F(v,10**6)) for v in (a,b,c,a+12,h+2)]
        hi = interpolate(double(A)+double(B)+double(C),powers,1)
        force = interpolate([(D,False),*double(B),*double(C),(H,False)],powers,1)
        lo = interpolate([*double(B),*double(C),(H,False)],powers[:-1],-1)
        targets = two["harmonics"](acb.zeta_zero(zero).imag,power,190,6)
        moments = arc["fifth"]["moment_boxes"](targets,powers)
        # Allow an additional error of 10^-12 in every measured harmonic.
        moments = [moments[0]] + [m+arb(0,q("1e-12")) for m in moments[1:]]
        assert value(hi,H)>0 and value(lo,D)>0
        assert mean(hi,moments)/value(hi,H) < q(F(1,75))
        assert mean(force,moments)>0
        assert mean(lo,moments)/value(lo,D) < q(F(1,75))
        for row,seeds in (
            (counts, (4874,8060,8111,8467,8927,9493)),
            ((1,1,21,8,48,1), (4828,7820,8086,8493,8929,9494)),
            ((1,1,21,13,58,1), (4641,6938,8098,8325,8940,9517)),
        ):
            weights = [q(F(n,sum(row))) for n in row]
            arc["upper_certificate"]([q(F(n,10000)) for n in seeds], targets,2**53,weights)
        print(zero,"excluded 75; exact uniform solutions at 76, 80 and 95.")


if __name__ == "__main__":
    verify()
