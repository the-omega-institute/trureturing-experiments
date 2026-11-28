"""验证两个含 theta23 八项误差的紧子水平上所需的统一常数。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb


def product3(values):
    return values[0]*values[1]*values[2]


def verify(replay):
    def load(name):
        return run_path(str(replay / f"ch39-{name}.py"))

    orbit = load("two-zero-orbit-condition")
    model, up, down = orbit["calibration"]()
    radius = arb(8)/5
    masses = [center+arb(0, (radius*width).upper()) for center, width in
              zip(model["ref"][:6], model["tol"][:6], strict=True)]
    for sector in (masses[:3], masses[3:]):
        assert 0 < sector[0] < sector[1] < sector[2]

    norm = orbit["frobenius"]
    rows = ((up, masses[2], ".000225", ".409", "2", "6700"),
            (down, masses[5], ".000437", ".407", ".034", "57"))
    for matrix, heaviest, sigma, operator, lower_scale, upper_scale in rows:
        sigma, operator, lower_scale, upper_scale = map(arb, (sigma, operator, lower_scale, upper_scale))
        assert 1/norm(matrix.inv()) > sigma and norm(matrix) < operator
        assert heaviest/operator > lower_scale
        assert arb(3).sqrt()*heaviest/sigma < upper_scale
    least_eigenvalue = (masses[3]/(57*arb(".407")))**2
    assert least_eigenvalue > arb("3.2e-13")

    tau = load("two-zero-sharp-quark-fit")["normalized_traces"](up,down)[0]
    kappa = (product3(masses[:3])/product3(masses[3:]))**(arb(2)/3)
    assert tau < 7200 and kappa < 52
    vtd_bound = (arb(7200)*52).sqrt()*masses[3]/masses[2]
    assert vtd_bound < arb(".01")
    angles = [v+arb(0,(w*radius).upper()) for v,w in zip(model["ref"][6:9],model["tol"][6:9],strict=True)]
    s12, s23 = angles[0].sin(), angles[1].sin()
    assert s12 < arb(".23") and s23 < arb(".044")
    assert angles[0].cos()*angles[1].cos() > arb(".97")
    free_s13 = (arb(".01")+arb(".23")*arb(".044"))/arb(".97")
    assert free_s13 < arb(".021") and angles[2].sin() < arb(".004")
    assert 1-arb(".021")**2 > arb(".9995")
    print("下型部门给出的最小归一化源本征值下界", least_eigenvalue.lower().str(12))
    print("相对谱直接迹上界与质量行列式因子上界", tau.upper().str(12),kappa.upper().str(12))
    print("留出第三角时 Vtd 模、第三角正弦的上界",vtd_bound.upper().str(12),free_s13.upper().str(12))
    load("shell-real-theta23-eight-upper")["verify"](replay)


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
