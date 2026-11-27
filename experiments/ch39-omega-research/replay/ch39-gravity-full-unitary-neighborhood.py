"""在十九维射影算子距离中，检验右侧两组荷六校准与轻上夸克窗口的距离。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    # 先前的四站下界在中心处容许正交补上任意上型酉矩阵。
    cube = run_path(str(root / "ch39-local-gravity-cube-root-flavor.py"))
    blocks = cube["verify"](root)
    parameters = ((".000020", ".000013"), (".000230", ".000040"))
    for (sites, center, box_radius), (radius_string, floor_string) in zip(blocks, parameters, strict=True):
        delta = arb(radius_string)
        floor = arb(floor_string)
        values = center.eig(algorithm="rump")
        inverse_gaps = []
        for i, value in enumerate(values):
            others = [z for k, z in enumerate(values) if i != k] + [1]
            gap = min(abs(value - z).lower() for z in others)
            assert gap > 2 * delta
            inverse_gaps.append(delta / (gap - delta))

        # 将各本征线及其共同正交补配准的极分解。平方和也界定补空间损失。
        squared = sum(x**2 for x in inverse_gaps)
        assert squared < 1
        alignment = (2 * squared).sqrt() + 1 - (1 - squared).sqrt()
        error = 2 * alignment + arb(4) / 3 * (delta / 2).asin()
        lower = floor.sqrt() - arb(".947555").sqrt() * error
        assert lower > arb(".00000629").sqrt()
        print("局部站、整空间射影半径、盒半径、轻上型振幅下界",
              sites, radius_string, box_radius.str(12), (lower**2).lower().str(12), flush=True)
    return blocks


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
