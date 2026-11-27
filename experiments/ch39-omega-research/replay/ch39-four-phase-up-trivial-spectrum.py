"""两条上型平凡线的谱区间已足以排除完整夸克中心，允许任改两条字符值。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    checker = run_path(str(root / "ch39-quark-trivial-eigenvalue-exclusion.py"))["verify"]
    # 第一盒向外舍入四相位根的两条平凡线。后三盒共用最轻上型目标作上界。
    boxes = (
        (".0000056259", ".0000058700", ".00400867", ".00401513"),
        ("0", ".00000615", ".005", ".008"),
        ("0", ".00000615", ".008", ".010"),
        ("0", ".00000615", ".010", ".018"),
    )
    for box in boxes:
        print("待排除的上型平凡块谱矩形", box, flush=True)
        checker(root, tuple(map(arb, box)))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
