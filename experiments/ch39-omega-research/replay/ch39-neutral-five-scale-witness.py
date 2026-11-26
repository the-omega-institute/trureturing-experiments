"""一次中性壳响应下接近尺度下界的七枚零点。"""
from pathlib import Path
from runpy import run_path

from flint import arb

replay = run_path(str(Path(__file__).with_name("ch39-seven-actual-neutral-times.py")))
indices = (8359, 143921293989490, 295873368963, 142270715263, 13219432687559, 295668865161, 142268163595)
scale = 2_903_616_618


def verify():
    gap = replay["verify_at"](scale, indices)
    assert arb("1.360899660e-8") < gap < arb("1.360899663e-8")


if __name__ == "__main__":
    verify()
