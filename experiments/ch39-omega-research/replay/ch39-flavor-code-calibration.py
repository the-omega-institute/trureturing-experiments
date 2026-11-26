"""Certify the finite input codes against the actual indexed zeta zeros."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb

root = Path(__file__).parent
finite = run_path(str(root/"ch39-rational-local-flavor.py"))
geo = run_path(str(root/"ch39-arc-relative-determinant.py"))


def verify():
    scale = finite["scale"]
    for (idx,codes),zero in zip(finite["records"](),acb.zeta_zeros(10815,1278),strict=True):
        phases = geo["phases"](zero.imag)
        for (sign,n),d in zip(codes,(724,804)):
            z = sign*phases[d]
            assert z.real > 0
            assert abs(scale*z.imag/(1+z.real)-n) < arb(1)/2
    assert arb("2.181") < arb("1.141").tan() < arb("2.182")
    assert arb("2.252") < arb("1.153").tan() < arb("2.253")
    print("certified both flavor clock codes for every zero 10815..12092")


if __name__ == "__main__":
    verify()
