"""在原十一壳插值处方上仅替换预登记的末端口支撑。"""
from pathlib import Path
from runpy import run_path


def verify(root, fifth, columns, levels):
    support = run_path(str(root / "ch39-terminal-port-support.py"))["verify"](root)
    experiment = run_path(str(root / "ch39-eleven-weak-eight-sites.py"))["experiment"]
    return experiment(root, fifth, columns, levels, sites=support,
                      pivot_bounds=(".02733", ".01643"), face_bound="1.4342e-8")


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    original = run_path(str(root / "ch39-fifth-full-positive-joint.py"))["verify"]
    verify(root, *original(root, return_operator="spectral"))
