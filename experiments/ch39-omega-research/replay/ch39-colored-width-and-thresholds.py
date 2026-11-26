"""轻色伙伴的直接 Z 宽度、低能强耦合与 W 质量层的电荷重分配。"""
from fractions import Fraction
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


root = Path(__file__).parent
broken = run_path(str(root / "ch39-broken-electroweak-vacuum.py"))
color = run_path(str(root / "ch39-all-light-color-assignments.py"))


def sm_width(alpha):
    # Dubovyk 等，式 (13) 和表 5。其余电弱输入取参考值，0.012 是拟合公式误差。
    lh = (arb("125.25") / arb("125.7")).log()
    dt = (arb("172.5") / arb("173.2"))**2 - 1
    ds = alpha / arb(".1184") - 1
    return (arb("2494.74") - arb("2.3")*lh + arb("19.9")*dt + arb("58.61")*ds
            - 4*ds**2 + 8*ds*dt + arb(0, ".012"))


def check_reallocation(spectrum, mz):
    # 每层实际标量实维数固定，外部规范场取零，牛顿权只依赖自旋与质量。
    for name in ("u", "d", "s", "c", "b", "W"):
        g, _, mass2 = spectrum[name]
        scalars, dirac = broken["minimal_partners"](g, 0)
        if name == "W":
            assert (scalars, dirac) == (2, 2) and mz/2 < mass2.sqrt() < mz
        else:
            assert (scalars, dirac) == (8 + 4, 0) and mass2.sqrt() < 10
        assert g + scalars - 4*dirac == 0
    # GUT 超荷归一。W 层一复单态给 dY²/5，实色伴随给 T(adj)/6。
    b1 = Fraction(1, 5)
    b3 = 5 * Fraction(3, 6)
    assert (15*b1, 6*b3) == (3, 15)
    color["check_tangencies"]([("8",)*5])
    return arb(b3.numerator)/b3.denominator


def verify():
    ctx.prec = 192
    v, mz, sin2 = arb("246.22"), arb("91.1876"), arb(".23129")
    spectrum = {name: (g, weight, mass2)
                for name, g, weight, mass2 in broken["light_spectrum"](arb(1))}
    masses = [spectrum[name][2].sqrt() for name in ("u", "d", "s", "c", "b")]
    width_unit = 1000*mz**3/(12*arb.pi()*v**2)*(sin2/3)**2

    def phase(mass):
        return (1 - 4*mass**2/mz**2)**(arb(3)/2)

    widths = [width_unit*(3*phase(a) + 6*phase(b)) for a, b in permutations(masses, 2)]
    width = min(w.lower() for w in widths).union(max(w.upper() for w in widths))
    assert len(widths) == 20 and arb("17.5920") < width < arb("17.7478")
    print("二十种分配的新增 Z 宽度/MeV", width.lower().str(10), width.upper().str(10))

    b3_shift = check_reallocation(spectrum, mz)
    log_interval = (mz/10).log()/(2*arb.pi())
    shift = b3_shift*log_interval
    assert arb(".87946") < shift < arb(".87947")
    print("固定 M_Z 边界，10 GeV 逆强耦合的正偏移", shift.str(11))
    for number, alpha in [(22027, arb(".1188549474", "3e-10")),
                          (28325, arb(".1187872682", "3e-10"))]:
        residual = sm_width(alpha) + width - arb("2495.2")
        assert arb("17.15") < residual < arb("17.38")
        print(number, "参考总宽加直接产生相对实测中心的差/MeV",
              residual.lower().str(9), residual.upper().str(9))
        inverse_sm = 1/alpha - arb(23)/3 * log_interval
        alpha_sm = 1/inverse_sm
        alpha_colored = 1/(inverse_sm + shift)
        assert arb(".174") < alpha_sm < arb(".176")
        assert arb(".151") < alpha_colored < arb(".153")
        print(number, "10 GeV 一圈标准模型及两种有色方案的强耦合",
              alpha_sm.str(8), alpha_colored.str(8))


if __name__ == "__main__":
    verify()
