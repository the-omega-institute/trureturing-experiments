"""第五壳六个大权层均低时，两个引力矩迫使第三及第四层同时处于高截止区。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    def necessary_region(fourth, third, trace):
        return fourth**2 > arb(".10") and third**2 > arb(".025") and trace-fourth-third < arb(".66")

    moments = run_path(str(root/"ch39-fourth-third-pivot-strong-gravity-cutoff.py"))
    return moments["verify"](root, number=2_306_888_100_090, certified_region=necessary_region)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
