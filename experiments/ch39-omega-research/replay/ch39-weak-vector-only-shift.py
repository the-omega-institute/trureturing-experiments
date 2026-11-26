"""两个中性读数的线性组合消去全部重门户，只读到带电弱规范场移质。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    ctx.prec = 256
    root = Path(__file__).resolve().parent
    vector = run_path(str(root / "ch39-em-colour.py"))["vector_log_moment"]
    mixed = run_path(str(root / "ch39-mixed-neutral.py"))["mixed_log_moment"]

    def logarithm(ratio):
        return 2*mixed(ratio)-7*vector(ratio)

    v, planck, mz, g2 = map(arb, ("246.22", "1.22089e19", "91.1876", ".63811"))
    mass = (g2*v/(2*planck))**2
    p0 = (mz/planck)**2
    min_p = (arb("1e13")/planck)**2
    x_min, x_max = arb("4e-15"), arb(1)
    factor = 1/(4*arb.pi())
    area = arb(22)/3
    large = factor*logarithm(mass/p0)
    tail = factor*logarithm(mass/min_p)
    delta = factor*area*mass/x_min
    lower_tv = 2*(large-tail)/delta-1
    error = tail+13*delta
    assert lower_tv > arb("3.56923074e20") and error < arb("8.17e-20")
    assert arb("1.07795564493676390410896") < large < arb("1.07795564493676390410897")
    print("轻规范及鬼的纯阈值、全能标总变差≤25时的误差上界、消去位移必需总变差",
          (-large).str(27),error.upper().str(18),lower_tv.lower().str(18))
    for name, m, bound in (("全窗","1e13","3.878e20"),("六式能标盒","1.4513765e14","3.571762e20")):
        p = (arb(m)/planck)**2
        eta = x_min/p
        assert 0 < eta < arb(".5")
        end_loss = 16*eta+16*eta*((1-eta)/eta).log()
        fraction = (p0/4+mass)/(2*x_min)+end_loss/area
        gap = delta*(1-fraction)-factor*area*mass/x_max
        assert gap > 0
        upper_tv = 1+2*large/gap
        assert upper_tv < arb(bound)
        print(name,"归一且纯阈值为零的两端节点代价上界",upper_tv.upper().str(18))


if __name__ == "__main__":
    verify()
