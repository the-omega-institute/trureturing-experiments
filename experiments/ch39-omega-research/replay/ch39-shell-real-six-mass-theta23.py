"""两个实谱方向给出原制备盒上的六质量与第二混合角共同见证。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb_mat, arb


def verify(replay):
    orbit = run_path(str(replay / "ch39-two-zero-orbit-condition.py"))
    model, up, down = orbit["calibration"]()
    shell = run_path(str(replay / "ch39-shell-pencil-six-mass-window.py"))
    frame, h = shell["shell_frame"](model)
    adj = model["adj"]
    rows = (("-.783717455278122", ".000418305420631677", "-.0000388045108273901"),
            ("-.265521789308497", "-.00270125750743603", ".00000704129514570978"),
            ("-.561502564281588", ".000693513247742955", ".000050831745050594"))
    factor = acb_mat([[arb(x) for x in row] for row in rows])
    assert abs(factor.det()) > 0
    levels, basis = model["eig"](factor * adj(factor))
    assert 0 < levels[0] < levels[1] < levels[2]
    ratio = h[0,1].real / h[1,2].real
    for i in (1, 2):
        assert abs(basis[0,i]) > 0
        y, z = [(basis[j,i]/basis[0,i]).real for j in (1, 2)]
        a = y*(1-ratio*z)
        b = z*(ratio+z)-y*y
        t = ratio+z-ratio*y*y
        minor = (ratio*b-t)*(-t)-a*a
        assert abs(minor) > 0
        print("谱线与末二阶子式", i+1, minor.str(12))
    pair = []
    for prepared, scale in ((up, "24.9334657026892"), (down, "3.52740222930554")):
        yukawa = arb(scale) * prepared * frame * factor
        pair.append(model["eig"](adj(yukawa)*yukawa))
    values = [v.sqrt() for levels, _ in pair for v in levels]
    values.extend(model["mixing_angles"](pair[0][1], pair[1][1]))
    errors = [abs(v-c)/w for v,c,w in zip(values,model["ref"],model["tol"],strict=True)]
    matched = [*range(6), 7]
    worst = max(errors[j].upper() for j in matched)
    assert worst < arb(".001105")
    assert all(errors[j] > 1 for j in (6, 8))
    print("六质量及第二混合角最大误差，以参考半宽计", worst.str(12))
    print("三个角与带号相位", [v.str(10) for v in values[6:]])
    print("两个留出角的误差，以参考半宽计", [errors[j].str(8) for j in (6, 8)])


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
