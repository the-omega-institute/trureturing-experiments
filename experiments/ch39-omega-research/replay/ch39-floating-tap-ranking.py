"""Certify the four-tap rankings on the floating-scale shortlist at the actual phases."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    ctx.prec = 160
    root_two = arb(2).sqrt()
    # Same 64-degree normalization and 45-bit absolute readout as ch39-long-born-return.
    cutoff = arb(".070")
    amplified = (32 * (arb(".0795438") / cutoff).acosh()).cosh()
    normalization = (32 * (arb(".079546") / cutoff).acosh()).cosh()
    leakage = 1 / (arb(".001") * amplified)**2
    port_floor = (amplified / normalization)**2 * arb(".001")**2 * (1 - leakage)
    rounding = 2 / (2**45 * port_floor)
    retained = ((1 - leakage) * (1 - rounding) / ((1 + leakage) * (1 + rounding)))**2
    coefficient_error = 4 * (1 - retained) * (arb(13) / 4)**2
    assert coefficient_error < arb(1) / 65536
    modulus_error = (3 + 3 * root_two) * coefficient_error

    weights = []
    dyadic = run_path(str(Path(__file__).with_name("ch39-dyadic-shortlist.py")))
    for index, g_bin, u_bin in dyadic["CELLS"]:
        zero = acb.zeta_zero(index)
        assert zero.real == arb(1) / 2
        height = zero.imag
        phase = root_two * acb(0, height * arb(2).log()).exp()
        scale = dyadic["SCALE"]
        assert (height * scale).floor().unique_fmpz() == g_bin
        assert (phase.real * scale).floor().unique_fmpz() == u_bin
        ideal = abs((phase - 1) * (phase - 2) * (9 + 4 * phase) / 13)
        denominator = height**2 + arb(1) / 4
        weights.append(((ideal - modulus_error) / denominator,
                        (ideal + modulus_error) / denominator))
    gaps = [lo - hi for (lo, _), (_, hi) in zip(weights, weights[1:])]
    assert all(gap > arb("4e-12") for gap in gaps)
    print("Tap error", coefficient_error.upper().str(8), "minimum certified weight gap",
          min(gaps).lower().str(8))
    print("Strict four-tap order: 10854, 10875, 10843, 10871, 10908, 10918")


if __name__ == "__main__":
    verify()
