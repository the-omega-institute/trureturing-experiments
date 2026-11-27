"""横跨 804 相位整圆，排除共享这只钟时的完整夸克窗。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify():
    certificate = run_path(str(Path(__file__).with_name("ch39-shared-724-quark-gaps.py")))
    total = 0
    for sign in (1, -1):
        for lo, hi, outer in ((0, ".515", True), (".515", ".95", False), (".95", arb.pi(), True)):
            tests = ((-4, 1, "2.5"), (2, 0, "14.7")) if outer else ((-4, 0, "11"), (2, 1, "4.15"))
            for charge, index, limit in tests:
                total += certificate["verify"](charge, sign, lo, hi, index, limit, shared_clock=804)
    assert arb("2.5")*arb("14.7")*(arb(9)/7)**2 < arb("61.32")
    assert arb("11")*arb("4.15")*(arb(9)/7)**2 < arb("77.12")
    certificate["cert"]["target_gaps"]()
    print("804 共同钟全域终格数", total)


if __name__ == "__main__":
    verify()
