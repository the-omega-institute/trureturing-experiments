"""从有向平凡轴、夸克谱间距与 CKM 导出四站物理超平面的投影距离下界。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    ctx.prec = 256
    standard = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))["standard_matrix"]
    angles = [arb(t).sin() for t in (".22704", ".04275", ".00377")]
    mix = standard(*angles, arb("1.139"))
    s = (1 - max(abs(z).upper()**2 for row in mix for z in row)).sqrt()
    up = list(map(arb, (".00000615", ".00311", ".8616")))
    down = list(map(arb, (".0000135", ".000268", ".01401")))
    gaps = [min(b - a for a, b in zip(m, m[1:])) for m in (up, down)]
    normal = arb(2).sqrt()*s/(1/gaps[0] + 1/gaps[1])
    print("不同荷任一本征线间距与平凡投影长度下界", s.str(14), normal.str(14))
    # 三角形证书的严格外向取整值。最后一对的两条例外满足 |ΔHx| > .02。
    for pair, chord in (("N−2/N−1", ".00189564"), ("N−1/N", ".30042033"),
                        ("N−2/N", ".01583171")):
        bound = normal * arb(chord)/2
        assert bound > 0
        print("壳对与两个实际三代投影的算子范数距离下界", pair, bound.str(12))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
