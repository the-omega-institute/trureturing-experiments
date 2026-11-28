"""在六质量和含 theta23 的两个八项问题上给出可复算的同时误差上界。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb_mat, arb


def verify(replay):
    model, up, down = run_path(str(replay / "ch39-two-zero-orbit-condition.py"))["calibration"]()
    frame, _ = run_path(str(replay / "ch39-shell-pencil-six-mass-window.py"))["shell_frame"](model)
    adj = model["adj"]
    cases = (
        (6, ((".8023655053776657", ".0007796139201564905", ".00001909686943567969"),
             (".2151543733143181", "-.003904900017689629", ".000001055060461370912"),
             (".5567029651652365", ".0003855215663375419", "-.00002793171788415042")),
         "27.00510257323013", "3.893566599733536", "1.539"),
        (8, (("-.8027178914011188", ".0007774645306128613", "-.00002022385393079161"),
             ("-.2133185741741397", "-.003927824232147434", "-.000001112292715808571"),
             ("-.5569014030659775", ".0003838977150687484", ".00002957672936151398")),
         "27.02700431743638", "3.821641857496311", "1.08"),
    )
    for extra, rows, up_scale, down_scale, bound in cases:
        factor = acb_mat([[arb(v) for v in row] for row in rows])
        assert abs(factor.det()) > 0
        pair = []
        for matrix, scale in ((up, up_scale), (down, down_scale)):
            yukawa = arb(scale)*matrix*frame*factor
            pair.append(model["eig"](adj(yukawa)*yukawa))
        readings = [lam.sqrt() for levels, _ in pair for lam in levels]
        readings.extend(model["mixing_angles"](pair[0][1], pair[1][1]))
        indices = [*range(6), 7, extra]
        deviations = [abs(readings[i]-model["ref"][i])/model["tol"][i] for i in indices]
        worst = max(q.upper() for q in deviations)
        print("第二个角的索引、八项最大偏差的上界", extra, worst.str(13))
        assert worst < arb(bound)


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
