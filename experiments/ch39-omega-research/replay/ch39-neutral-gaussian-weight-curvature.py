"""真实辐边的中性高斯行列式在对数高度上的曲率界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

shell = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))


def verify():
    ctx.prec = 256
    edges = shell["geo"]["edges"]
    pole, tuned = shell["soft"]["pole"], shell["soft"]["star"]["orbit"]
    offset = shell["soft"]["point"]["study"]["offset"]
    assert offset > 1
    size = 19
    fixed, adjusted = arb_mat(size, size), arb_mat(size, size)
    for i, j, _, _ in edges:
        for lap in [fixed] + ([adjusted] if pole in (i, j) and (i in tuned or j in tuned) else []):
            lap[i,i] += 1
            lap[j,j] += 1
            lap[i,j] -= 1
            lap[j,i] -= 1
    assert len(edges) == 135 and adjusted.trace() == 6
    assert max(fixed[i,i] for i in range(size)) <= 18
    d = arb(23)/(25*2**20)
    vertical = (d/offset+36*d*128**2)/(1-d/offset)**2
    assert 1-vertical > arb(".48249")
    first = acb.zeta_zero(8359).imag
    distortion = -3*(1-d/(first**2+offset)).log()
    minimum_power = arb(259).log()/distortion
    assert minimum_power > arb("1.51171e14")

    # 直接求三个有限参数例的二阶导数，交叉核验非交换项与正号。
    eye = arb_mat([[int(i == j) for j in range(size)] for i in range(size)])
    for b, mass, t in [(arb(1)/128, 0, arb(1)), (arb(1), 3, first**2),
                        (arb(2903600000), 0, arb(10)**26)]:
        x = b*b/t
        delta = d/(t+offset)
        u = t/(t+offset)
        slope = -delta*u
        curvature = delta*u*(2*u-1)
        k = fixed-adjusted*delta+eye*mass
        inv = (k+eye*x).inv()
        dt = -eye*x-adjusted*slope
        ddt = eye*x-adjusted*curvature
        second = (inv*ddt-inv*dt*inv*dt).trace()
        base = (inv*inv*k).trace()*x
        assert base > 0 and second > arb(".48249")*base
    print("相对固定谱曲率的全域下界", (1-vertical).str(13))
    print("允许任意小核尺度时，由三条辐边提供所需校正的总逆幂下界", minimum_power.str(16))


if __name__ == "__main__":
    verify()
