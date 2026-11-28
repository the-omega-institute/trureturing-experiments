"""校准四源中，附加项限于有向交换子支撑平面时的最小正迹份额。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def verify():
    replay = Path(__file__).resolve().parent
    def load(name):
        return run_path(str(replay/("ch39-"+name+".py")))

    orbit = load("two-zero-orbit-condition")
    model,_,_ = orbit["calibration"]()
    adj = model["adj"]
    neighbour = load("shell-commutator-next-zeros")
    at_zero = load("quadratic-gram-source")["at_zero"]
    up,down = [at_zero(model,q,n) for q,n in ((-4,2345666),(2,964890))]
    parts,_ = orbit["orbit"](model,up,down)
    chart,_ = load("shell-pencil-six-mass-window")["shell_frame"](model)
    axes,_ = neighbour["shell_axes"](model)
    words = load("two-zero-quadratic-source")
    equations = words["orbit_coefficients"](parts,adj,axes)
    trig = words["trigonometric_columns"]
    centers,_ = neighbour["certify_centers"](equations,trig,neighbour["CENTERS"][0])
    refine = load("four-quadratic-source-neighbours")["refine"]
    # 依第一相位递增排序，对应四个精确的中央十项见证。
    windows = (("6.16","6.18"),("7.61","7.63"),("49.34","49.37"),("13.30","13.33"))
    for (x0,y0),(lo,hi) in zip(centers,windows,strict=True):
        x,y = refine((x0,y0,0),arb(".02"),equations,trig)
        matrix = adj(chart)*orbit["phase_sum"](parts,orbit["phases_at"](x,y))
        covariance = matrix*adj(matrix)
        real = (covariance+covariance.conjugate())/2
        amplitude = abs(covariance[0,1].imag)
        schur = acb_mat([[real[i,j]-real[i,2]*real[2,j]/real[2,2] for j in (0,1)]
                         for i in (0,1)])
        levels,_ = model["eig"](schur)
        small,large = levels
        assert 0 < small < amplitude and small*large>amplitude**2
        share = (small+amplitude**2/small)/covariance.trace().real
        assert arb(lo) < 10**6*share < arb(hi)
        print("第一相位与最小有向正项的百万分迹份额",x0.str(7),(10**6*share).str(9))


if __name__=="__main__":
    verify()
