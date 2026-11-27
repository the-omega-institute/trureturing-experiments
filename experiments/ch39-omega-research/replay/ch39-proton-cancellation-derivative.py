"""第 177 枚正零点处，八站固定非负电流相消时的高度导数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def verify():
    model = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))
    ctx.prec = 192
    geo, soft = model["geo"], model["soft"]
    sites = (0, 2, 3, 5, 6, 9, 11, 13)
    targets = (3, 5)
    height = acb.zeta_zero(177).imag
    laplacian = model["laplacian"](5, height)
    green = laplacian.inv()

    speeds = {724: 5*geo["a"], 804: 5*geo["b"]}
    derivative = acb_mat(19,19)
    for i,j,distance,_ in geo["edges"]:
        if distance in speeds:
            assert soft["pole"] not in (i,j)
            derivative[i,j] = acb(0,speeds[distance])*laplacian[i,j]
            derivative[j,i] = derivative[i,j].conjugate()
    frozen_slope = -green*derivative*green

    denominator = height**2 + soft["point"]["study"]["offset"]
    spoke_slope = 2*height*(arb(23)/25)*arb(2)**-20/denominator**2
    pole = soft["pole"]
    for j in soft["star"]["orbit"]:
        derivative[pole,pole] += spoke_slope
        derivative[j,j] += spoke_slope
        derivative[pole,j] -= spoke_slope
        derivative[j,pole] -= spoke_slope
    full_slope = -green*derivative*green

    forward = [(6,j) for j in targets]
    backward = [(j,6) for j in targets]
    other = [(i,j) for i in sites for j in sites
             if (i,j) not in forward and (i,j) not in backward]
    assert len(other) == 60
    assert all(green[i,j].real > arb("-.00075") for i,j in forward+backward)
    assert all(green[i,j].real > arb(".006") for i,j in other)
    for slope in (frozen_slope,full_slope):
        assert all(arb(".0195") < slope[i,j].imag < arb(".0204") for i,j in forward)
        assert all(arb("-.0204") < slope[i,j].imag < arb("-.0195") for i,j in backward)
        assert all(abs(slope[i,j].imag) < arb(".0065") for i,j in other)

    special_mass = arb(8)/9
    reverse_mass = (1-special_mass)**2/(2*special_mass)
    lower = arb(".0195")*special_mass-arb(".0399")*reverse_mass-arb(".0065")*(1-special_mass)
    assert lower > arb(".0163")

    # 取前一实验中精确相消的三个实权，独立核对导数的符号和量级。
    x,y,z = [green[6,j] for j in (3,5,13)]
    areas = [(u.conjugate()*v).imag for u,v in ((y,z),(z,x),(x,y))]
    assert all(w < 0 for w in areas)
    weights = [w/sum(areas) for w in areas]
    slope = sum(w*full_slope[6,j].imag for w,j in zip(weights,(3,5,13),strict=True))
    assert arb(".01922") < slope < arb(".01923")
    print("统一导数下界、面积权处的虚部导数", lower.str(12), slope.str(12))


if __name__ == "__main__":
    verify()
