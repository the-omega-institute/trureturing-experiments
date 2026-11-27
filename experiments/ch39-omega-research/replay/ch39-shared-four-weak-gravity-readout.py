"""按事前登记在同长四门的两个三弱解盒中逐项读取引力。"""
from pathlib import Path
from runpy import run_path

def verify(root):
    src = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    boxes = run_path(str(root/"ch39-shared-four-gates-weak-boxes.py"))
    readings = boxes["experiment"](root,*src["verify"](root,return_operator="spectral"),read_gravity=True)
    assert len(readings)==6
    for k,w in readings:
        assert k < -32 and w > 14
    for first,second in zip(readings[:3],readings[3:],strict=True):
        assert first[0] < second[0] and first[1] > second[1]


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
