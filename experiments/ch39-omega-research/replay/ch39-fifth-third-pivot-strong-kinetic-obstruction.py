"""第五壳第三枢轴在原能标处有强权层超过阈值时的规范动力学排除。"""

from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    core = run_path(str(root / "ch39-fourth-third-pivot-global-kinetic-obstruction.py"))
    tables = core["build_grid"](root, (0,), number=2_306_888_100_090)[:5]
    results = core["verify_branches"](
        root,
        tables,
        "ch39-fifth-third-pivot-strong-kinetic-duals.txt",
        minimum_deficit=arb(".1305"),
    )
    assert [count for count, margin in results] == [87, 441, 90, 242, 38, 213]
    return results


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
