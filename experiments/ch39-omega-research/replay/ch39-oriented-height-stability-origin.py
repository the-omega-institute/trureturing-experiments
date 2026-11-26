"""Check the 46 actual-height columns against the rational reconstruction matrix."""
import json
from pathlib import Path
from runpy import run_path
from flint import arb


def main():
    root = Path(__file__).resolve().parent
    checker = run_path(str(root/"ch39-oriented-height-stability.py"))
    model = run_path(str(root/"ch39-five-tick-record-origin.py"))["model"]
    _, heights, weights = model()
    logs = [arb(x).log() for x in (2,3,arb(3)/2,6)]
    target = checker["matrix"](json.loads((root/"ch39-five-tick-record.json").read_text()))
    bound = arb(2)**-checker["BITS"]
    rows = [[w*g*g/(g*g+j*j)*(5*g*log).sin()*7/68 for g,w in zip(heights,weights)]
            for j in checker["OFFSETS"] for log in logs][:46]
    for row, exact in zip(rows,target):
        for val, fraction in zip(row,exact):
            assert abs(val-arb(fraction.numerator)/fraction.denominator) < bound
    print("All 2,116 coefficients from the 46 original height intervals agree within 2^-144.")


if __name__ == "__main__":
    main()

