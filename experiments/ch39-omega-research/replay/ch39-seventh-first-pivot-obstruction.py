"""第七壳第一枢轴的荷六五式在原匹配能标无解。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    endpoints = run_path(str(root/"ch39-fifth-arbitrary-pivot-obstructions.py"))
    first, = endpoints["verify"](root, number=2_306_888_100_092, pivots_to_check=(1,))
    assert first[0] > arb("366.37")
    return first[0]


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
