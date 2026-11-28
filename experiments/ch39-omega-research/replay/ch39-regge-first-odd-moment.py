"""Certify the tail inputs and the least Chebyshev witness degree for those inputs."""

from pathlib import Path
from runpy import run_path

from flint import arb, ctx, fmpq


def main():
    ctx.prec = 160
    support = run_path(str(Path(__file__).with_name("ch39-regge-shape-support.py")))
    width = fmpq(3, 50000)
    action, volume, fatness = support["chain"](1, width, "odd-moment")
    scale = arb(801) / 100
    total_volume = 6144 * volume
    assert fatness > fmpq(1, 16) and total_volume < scale**2
    b = arb(14402) / 100
    assert 6144 * action / total_volume.sqrt() < -b
    mass = (2 * arb(width) / scale)**3840 / support["enclosing_volume"]()
    exponent = 16314
    assert mass > arb(10)**(-exponent)
    a = 5 * arb(6144).sqrt() / 4
    threshold = (arb(10)**exponent - 1).acosh() / (b / a).acosh()
    assert 40176 < threshold < 40177
    print("Upper shape endpoint:", a)
    print("Negative tail exceeds 10^-16314 below -144.02.")
    print("Chebyshev threshold:", threshold)


if __name__ == "__main__":
    main()
