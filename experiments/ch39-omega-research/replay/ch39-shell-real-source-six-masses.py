"""用壳循环与极点的实结构复现六个参考质量，两个标尺和源均预先写在脚本中。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def verify():
    root = Path(__file__).resolve().parent
    model, up, down = run_path(str(root/"ch39-two-zero-orbit-condition.py"))["calibration"]()
    frame, _ = run_path(str(root/"ch39-shell-pencil-six-mass-window.py"))["shell_frame"](model)
    rows = ((-218953425, 1000000000000, -149653820),
            (80992301, 440865289757, 1133578394),
            (255353826, 717618895897, -487865564))
    factor = acb_mat([[arb(n)/10**12 for n in row] for row in rows])
    assert abs(factor.det()) > 0
    masses = []
    for prepared, scale in ((up, "17.1011054860"), (down,"2.4193378623")):
        masses.extend(model["singular"](arb(scale)*prepared*frame*factor))
    deviation = max((abs(v-m)/w).abs_upper() for v,m,w in
                    zip(masses, model["ref"][:6], model["tol"][:6], strict=True))
    assert deviation < arb(".0016")
    print("实源六质量最大偏差所占参考半宽", deviation.str(12))


if __name__ == "__main__":
    verify()
