"""A common rational portal occurs on a positive-density set of heights."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

finite = run_path(str(Path(__file__).with_name("ch39-shell-gw-finite-mass-eigenvalue.py")))
q = finite["q"]


def brackets():
    left, right = q("-.00675"), q("-.0035")
    band = q(125)**2 / (2*finite["unit"])
    m, _, error, _ = finite["face"](left, 5)
    left_upper = (m+error).upper()
    m, heavy, error, _ = finite["face"](right, 5)
    right_lower = (min(m.lower(), heavy)-error).lower()
    print("signed faces, in 1e24 GeV²", left_upper.str(8), right_lower.str(8))
    assert left_upper < -band and right_lower > band


def grid():
    left, right = Fraction(-675, 100000), Fraction(-35, 10000)
    step = Fraction(32, 10**23)
    spacing_count = (right-left)/step
    assert spacing_count.denominator == 1
    points = int(spacing_count)+1
    assert step < Fraction(125**2, 47990*10**21)
    # The zero-window's normalized torus area is 9e-16 / pi².
    density = q(9)/(10**16*points*q(1).pi()**2)
    assert density > q("8.97e-36")
    print("rational grid size", points, "relative density at least", density.lower().str(10))


if __name__ == "__main__":
    brackets()
    grid()
