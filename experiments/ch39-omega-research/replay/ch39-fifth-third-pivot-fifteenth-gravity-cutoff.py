"""第五壳第三枢轴在六个强权层均低时必有第十五层平方截止超过 0.001。"""

from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    gravity = run_path(str(root / "ch39-fourth-third-pivot-strong-gravity-cutoff.py"))
    return gravity["verify"](
        root, number=2_306_888_100_090, capped_layer=(14, arb(".001"))
    )


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
