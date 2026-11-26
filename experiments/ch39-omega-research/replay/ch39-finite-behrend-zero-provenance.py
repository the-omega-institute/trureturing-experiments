"""Independent calibration of the finite input table against the five specified zeta zeros."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

record = run_path(str(Path(__file__).with_name("ch39-finite-behrend-zero-record.py")))
ctx.prec = 160


def ball(v):
    return arb(v.numerator)/v.denominator


def main():
    cols,weights = record["basis_and_weights"]()
    residual = arb(0)
    for row, col in zip(record["records"],cols):
        h = acb.zeta_zero(row[0]).imag
        for base,stored in zip((2,3),col[1:3]):
            true = acb(0,h*arb(base).log()).exp()
            residual = max(residual,abs(true-acb(*(ball(v) for v in stored))).upper())
    assert residual < 1/arb(2)**22
    # The product rotor error is bounded by twice the base error.
    print("Maximum chord error for the ten primitive records",residual.str(10))


if __name__ == "__main__":
    main()
