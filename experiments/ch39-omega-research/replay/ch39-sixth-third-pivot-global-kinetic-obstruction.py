"""第六壳实际第三枢轴在原能标处的五式与全程正动能排除。"""

from pathlib import Path
from runpy import run_path

from flint import arb


def verify_branches(root, *, portal=None):
    core = run_path(str(root / "ch39-fourth-third-pivot-global-kinetic-obstruction.py"))
    # 源谱在 build_grid 中固定为 N+1。参数 number 选取目标壳的实际谱权和破缺背景。
    tables = core["build_grid"](root, (0,), number=2_306_888_100_091, portal=portal)[:5]
    return core["verify_branches"](
        root,
        tables,
        "ch39-sixth-third-pivot-global-kinetic-duals.txt",
        minimum_deficit=arb(".1968"),
    )


def verify(root, *, portal=None):
    gravity = run_path(str(root / "ch39-sixth-third-pivot-gravity-cutoff.py"))
    assert gravity["verify"](root, portal=portal) == (1, 10, 16)
    return verify_branches(root, portal=portal)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
