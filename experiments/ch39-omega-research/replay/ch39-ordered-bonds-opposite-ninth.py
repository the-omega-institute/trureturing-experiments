"""固定同一门序与第七、第八壳四个引力端点，认证第九壳的另一组引力符号。"""
from pathlib import Path
from runpy import run_path

from flint import arb


PARAMETERS = ("-1.921875", "1.15625", "-2.79795713653470", "-.359375",
              "-.203125", "-.109375", "1.15625", "-2.50710126346897",
              "-1.953125", ".703125", ".71497631174600", "-3.00899521637767")
ACTIVE = (2, 7, 10, 11)


def experiment(root, fifth, columns, levels):
    ordered = run_path(str(root / "ch39-two-shell-ordered-bond-gates.py"))
    box, _ = ordered["experiment"](root, fifth, columns, levels, parameters=PARAMETERS,
                                   active=ACTIVE, half_width=".00005", maximum_ratio=".769")
    held = run_path(str(root / "ch39-ordered-bond-ninth-holdout.py"))
    k, w = held["read_ninth"](root, fifth, columns, levels, box, parameters=PARAMETERS, active=ACTIVE)
    assert arb("-1.002893") < k < arb("-1.001641")
    assert arb(".003965") < w < arb(".004002")
    return k, w


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))
