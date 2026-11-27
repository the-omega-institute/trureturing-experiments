"""第五壳任意相干正解经任意强度的第六壳非选择测量后仍排除第七壳真空。"""

from pathlib import Path
from runpy import run_path


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    verify = run_path(str(root / "ch39-fifth-seventh-common-positive-operator.py"))[
        "verify"
    ]
    verify(root, sixth_measurement=True)
