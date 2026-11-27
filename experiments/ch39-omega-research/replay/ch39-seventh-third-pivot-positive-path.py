"""第七壳原生第三枢轴的荷六五项端点与全程正动能。"""
from pathlib import Path
from runpy import run_path

from flint import arb


CENTERS = (
    "1.6786997691972188e-34", ".21107063869554935", "2.0893539188522954e-26",
    "1.9251881555446e-26", ".05261846065429505", "2.297570158689684e-21",
    ".2460623899635208", "6.03804636069189e-25", "1.5229371614203173e-23",
    "5.236611405688443e-8", "3.5729317177099885e-24", "1.0482125343242282e-24",
    "2.614868378930508e-25", "1.2377133898300254e-23", "1.6484213228062992e-12",
    "6.700963689093833e-34", "3.94193610796293e-34", "3.925628567131245e-34",
    "3.936732812420191e-34",
)


def verify(root):
    solve = run_path(str(root / "ch39-fifth-third-pivot-em-colour-path.py"))["solve"]
    data, species, weights, centers, boxes, weak, lib = solve(
        root, CENTERS, (17, 0, 15, 1, 4, 6), number=2_306_888_100_092
    )
    threshold = arb("1e-6")
    high = (1, 4, 6)
    assert all(x > threshold if j in high else x < threshold for j, x in enumerate(boxes))
    path = lib("fourth-third-pivot-neutral-positive-path")["path"]
    assert weak < arb("7.600947")
    leaves, floors = path(data, species, weights, centers, boxes, weak, lib)
    assert floors[2] > arb("2.31401") and floors[3] > arb("4.71292")
    return leaves, floors


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
