"""固定右内壳源谱，对所有十九维酉输运求两个非平凡字符振幅的精确极小值。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def shell_data(root):
    training = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    spectral = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    source_number = training["CENTRAL"] + 1
    _, labels, _ = training["verify_nodes"](
        root, ((source_number, ".0294"),), ".1", (), ((6, 7, 10, 13), (6, 7, 10, 12))
    )
    ctx.prec = 384
    _, target = spectral["spectrum"](shell, -4, source_number + 1)
    assert (shell["soft"]["pole"], *shell["boundary_sites"]()) == (14, 8, 7, 6)
    weights = labels[source_number, -4]
    assert len(weights) == 19 and min(weights) > 0
    return weights, target


def ordered_product(weights, probabilities):
    # 上下端分别排序可处理重复的源权与恰为零的字符投影。
    low = sorted(w.lower() for w in weights)
    high = sorted(w.upper() for w in weights)
    p_low = sorted((max(arb(0), p.lower()) for p in probabilities), reverse=True)
    p_high = sorted((p.upper() for p in probabilities), reverse=True)
    lower = sum(x * y for x, y in zip(low, p_low, strict=True))
    upper = sum(x * y for x, y in zip(high, p_high, strict=True))
    return lower.union(upper)


def verify(root):
    weights, target = shell_data(root)
    omega = acb(-1, arb(3).sqrt()) / 2
    limits = (arb(".00000629500131095"), arb(".0000154660569880"))
    for k, limit in enumerate(limits, start=1):
        probabilities = [
            abs(target[8, j] + omega**k * target[7, j] + omega**(2*k) * target[6, j])**2 / 3
            for j in range(19)
        ]
        assert sum(probabilities).contains(1)
        assert sum(p > 0 for p in probabilities) == 6
        minimum = ordered_product(weights, probabilities)
        assert minimum > limit > arb(".00000629")
        print("字符指数、任意十九维酉输运上的最小振幅", k, minimum.str(16))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
