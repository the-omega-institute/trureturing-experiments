"""保持已登记换基与标尺，在两个实际零点的相邻编号读出六质量与混合。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def verify():
    calibration = run_path(str(Path(__file__).with_name("ch39-two-zero-quark-window.py")))
    graph = calibration["shell"]
    source = calibration["sources"]
    adjoint = calibration["adj"]
    eig = calibration["eig"]
    reference = calibration["ref"]
    widths = calibration["tol"]

    def classify(values, targets, tolerances):
        errors = [abs(y-t)/w for y,t,w in zip(values,targets,tolerances,strict=True)]
        assert all(z < 1 or z > 1 for z in errors)
        return sum(z < 1 for z in errors), [z.str(9) for z in errors]

    sectors = []
    for sector, (q, n, scale) in enumerate(((-4,2345666,"7.06539940386823251"),(2,964890,"1"))):
        results = []
        for m in range(n-1,n+2):
            height = acb.zeta_zero(m).imag
            laplacian = graph["laplacian"](q,height)
            rows = []
            for frame in calibration["frames"]:
                levels,basis = eig(adjoint(frame)*laplacian*frame)
                assert levels[0] < arb("8.685") and levels[1] > arb("9.6")
                first = acb_mat([[basis[i,0]] for i in range(basis.nrows())])
                projected = adjoint(frame*first)*source/(1+(levels[0]/9)**512)
                rows.append(projected.entries())
            y = acb_mat(rows)*calibration["K"]*arb(scale)
            squares,vectors = eig(adjoint(y)*y)
            masses = [x.sqrt() for x in squares]
            span = slice(3*sector,3*sector+3)
            count, errors = classify(masses,reference[span],widths[span])
            results.append((m,count,vectors))
            print("部门、编号、三质量、偏差的参考半宽数",q,m,[z.str(10) for z in masses],errors)
        sectors.append(results)
    for n,cu,u in sectors[0]:
        for m,cd,d in sectors[1]:
            angles = calibration["mixing_angles"](u,d)
            cv, _ = classify(angles,reference[6:],widths[6:])
            print("上下编号、六质量命中数、混合四项命中数",n,m,cu+cd,cv)


if __name__ == "__main__":
    verify()
