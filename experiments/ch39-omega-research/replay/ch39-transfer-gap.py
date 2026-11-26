"""Finite-slice positive transfer gap for unamplified zero weights."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
from flint import arb

ROOT = Path(__file__).resolve().parent
correlation = run_path(str(ROOT/"ch39-haar-correlation-length.py"))


def verify():
    moments = correlation["singlet_moments"]()
    high = correlation["influence"](1/(F(11134)**2+F(1,4)),moments)
    bounds = [(F(".3893"),F(".151555"),"1.334"),
              (F(".000004866"),F("2.368e-11"),"17.30"),
              (high,F("2.986e-13"),"20.39")]
    for rate, ratio, ea in bounds:
        value = arb(rate.numerator)/rate.denominator
        assert 0 < rate**2 < ratio and -arb(2).sqrt()*value.log() > arb(ea)
        print("Influence, two-layer spectral ratio, lower energy times spacing:",
              value.str(10),float(ratio),ea)


if __name__ == "__main__":
    verify()
