"""五壳十五式见证的预登记第十二壳读数。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def experiment(root, fifth, columns, levels):
    readout = run_path(str(root/"ch39-cyclic-nine-tenth.py"))["readout"]
    values = readout(root, fifth, columns, levels, "five-trivial-fifteen", 2_306_888_100_097)
    enclosures = (("1.24228", "1.24291"), ("-.040502", "-.040461"),
                  ("-.074256577194664", "-.074256577194663"))
    for value,(lo,hi) in zip(values,enclosures,strict=True):
        assert arb(lo)<value<arb(hi)
    return values


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    joint = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root, *joint["verify"](root, return_operator="spectral"))
