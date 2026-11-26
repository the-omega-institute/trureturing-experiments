"""Test gauge-coordinate and projective qubit volume on matched input segments."""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb


def uniform_g(x):
    # dg_i/dX_i in absolute value, with the common factors 1/2 omitted.
    return 1/prod(v*v.sqrt() for v in x)


def projective(x):
    # Fubini--Study line element dg/(1+g^2), with X=g^-2.
    return 1/prod(v.sqrt()*(1+v) for v in x)


if __name__ == "__main__":
    cert = run_path(str(Path(__file__).with_name("ch39-pointwise-volume-winner.py")))
    for prior in (uniform_g,projective):
        print(prior.__name__,flush=True)
        cert["verify"](arb(".000097"),prior)
