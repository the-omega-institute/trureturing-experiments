"""完整联合见证固定十九个物理截止及两型四站块后的局部维数。"""
from pathlib import Path
from runpy import run_path


def verify(root):
    joint = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    joint["verify"](root, certify_rank=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
