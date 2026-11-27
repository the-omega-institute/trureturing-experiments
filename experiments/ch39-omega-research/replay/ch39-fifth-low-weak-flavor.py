"""六个强权层均低时，完整夸克与第三枢轴的弱差零点可以同时成立。"""
from pathlib import Path
from runpy import run_path


if __name__ == "__main__":
    root=Path(__file__).resolve().parent
    run_path(str(root/"ch39-fifth-six-positive-flavor.py"))["verify"](root, weak_zero=True)
