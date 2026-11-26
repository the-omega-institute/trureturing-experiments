"""普朗克上限放开时，第 2、30 枚共用非负轮廓的最佳大截止误差系数。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify():
    root = Path(__file__).resolve().parent
    fields = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    _, _, rows = fields(root, (2, 30))
    first, second = rows
    assert sum(g for g, _, _ in first) == sum(g for g, _, _ in second) == -310
    assert sum(a for _, a, _ in first).overlaps(sum(a for _, a, _ in second))
    assert all(g == other_g for (g, _, _), (other_g, _, _) in zip(first, second, strict=True))
    difference = sum(g*(h-z) for (g, _, z), (_, _, h) in zip(first, second, strict=True))
    c = difference/310
    assert arb("1.6051e-8") < c < arb("1.6054e-8")
    # W'_30-W'_2 的正极限除以 W'_2/x 的正极限给出 c。
    print("两个谱的有符号质量矩增量，三个及任意多个截止的渐近系数",
          difference.str(15), c.str(15), (c/2).str(15))


if __name__ == "__main__":
    verify()
