"""第六壳三枚枢轴的正动能排除对十六个可各自选择的标量门户系数的区间重放。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    model = run_path(str(root/"ch39-sixth-third-pivot-global-kinetic-obstruction.py"))
    # 每次用区间平移相应的标量质量。各层在同一盒中可独立选值。
    portal = arb(0, 1)
    endpoints = run_path(str(root/"ch39-fifth-arbitrary-pivot-obstructions.py"))
    first, second = endpoints["verify"](root, number=2_306_888_100_091, portal=portal)
    assert first[0] > 429 and second[0] > 237 and second[1] > 3699
    return model["verify"](root, portal=portal)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
