"""第七壳实际第二枢轴满足三个规范端点时，在三个固定能标各有负色动能。"""

from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    core = run_path(str(root / "ch39-fourth-third-pivot-global-kinetic-obstruction.py"))
    tables = core["build_grid"](root, (0,), number=2_306_888_100_092, pivot=2)[:5]
    lower, _, high, weights, _ = tables
    for energy_index, exponent, deficit, count in (
        (1, 5, "14.6877", 17),
        (2, 7, "15.5983", 17),
        (3, 9, "2.4812", 19),
    ):
        filename = f"ch39-seventh-second-pivot-colour-{exponent}-duals.txt"
        # 解码后三个端点在先，此后交替为电磁与色。每棵树仅使用指定能标的色读数。
        allowed = 3 + 2 * energy_index + 1
        for _, _, direction in core["source_boxes"](
            root, high, len(lower)-1, weights, filename=filename, whole_domain=True
        ):
            assert all(v == 0 for i, v in enumerate(direction) if i >= 3 and i != allowed)
        results = core["verify_branches"](
            root, tables, filename, minimum_deficit=arb(deficit), whole_domain=True
        )
        assert len(results) == 1 and results[0][0] == count
        print("检验能标/GeV、色逆耦合严格上界", 5.6*10**exponent, "-"+deficit, flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
