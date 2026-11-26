"""Push quaternion Haar preparations to the three low-energy gauge channels."""
from fractions import Fraction
from math import comb, prod
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

cert = run_path(str(Path(__file__).with_name("ch39-pointwise-volume-winner.py")))


def radial_rotors(x):
    # Three independent rotors, g = ||Im q||/|Re q| for each.
    return 1/prod(v.sqrt()*(1+v)**2 for v in x)


def common_rotor(x):
    # One isotropic rotor, g_i = |q_i/q_0|, with a shared denominator.
    return 1/(prod(v*v.sqrt() for v in x)*(1+sum(1/v for v in x))**2)


def verify():
    ctx.prec = 112
    # Check the first 64 nonconstant moments against the rational scalar-axis recurrence.
    moment = Fraction(1)
    for n in range(64):
        moment *= Fraction(2*n+1,2*(n+2))
        assert moment == Fraction(comb(2*n+2,n+1),(n+2)*4**(n+1))
    v = cert["volume"]
    center = v["gauge"]["centers"]
    g = [c*v["q"](".999").union(v["q"]("1.001")) for c in center]
    box_size = prod(v["q"](".002")*c for c in center)
    radial_mass = box_size*prod(4*a**2/(arb.pi()*(1+a**2)**2) for a in g)
    shared_mass = box_size*8/(arb.pi()**2*(1+sum(a**2 for a in g))**2)
    assert v["q"]("3.89e-11") < radial_mass < v["q"]("3.98e-11")
    assert v["q"]("2.70e-10") < shared_mass < v["q"]("2.73e-10")
    assert 1/(1+sum(a**2 for a in g)) > v["q"](".36")
    print("Unconditioned probability of the low-energy gauge box, independent/common",
          radial_mass.lower().str(7),radial_mass.upper().str(7),
          shared_mass.lower().str(7),shared_mass.upper().str(7),flush=True)
    for density in (radial_rotors,common_rotor):
        print(density.__name__, flush=True)
        cert["verify"](arb(".000097"),density)


if __name__ == "__main__":
    verify()
