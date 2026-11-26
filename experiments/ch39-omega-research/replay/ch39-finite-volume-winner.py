"""Certify the volume-weighted maximum at every residual half-width up to 2^-21."""
from pathlib import Path
from runpy import run_path

from flint import arb


if __name__ == "__main__":
    volume = run_path(str(Path(__file__).with_name("ch39-volume-contenders.py")))
    volume["verify"](arb(2)**-21)
