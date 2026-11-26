"""在整个联合相位窗，把十二面质子因子同时压低的极限与实际零点密度。"""
from pathlib import Path
from runpy import run_path

from flint import arb

if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    density = run_path(str(root / "ch39-joint-zero-vacuum-flavor-density.py"))["verify"]()[-2]
    run_path(str(root / "ch39-balanced-proton-residual.py"))["verify"](joint_window=True)
    lower = density*arb(".293")
    assert lower > arb("5.3436e-14")
    print("实际零点通过原作用与味窗且新因子超阈值的相对下密度",lower.lower().str(14))
