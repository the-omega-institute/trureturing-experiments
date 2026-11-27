"""第七壳实际第三枢轴的四个正权层中至少一层有高截止。"""

from pathlib import Path
from runpy import run_path


def verify(root):
    gravity = run_path(str(root / "ch39-sixth-third-pivot-gravity-cutoff.py"))
    return gravity["verify"](
        root, number=2_306_888_100_092, positive_layers=(1, 7, 10, 16), maximum=".896034"
    )


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
