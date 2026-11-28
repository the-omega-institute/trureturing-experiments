"""第 6636 枚作为两荷共用零点时，纯类内等权可实现其校准多项式。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


# 只用类顶点。五个整数依次为上型轨道数、极点标志、下型轨道数、
# 极点标志、共同三轨道数，末列给出上型与下型的边界。
CLASSES = (
    (1,1,1,1,1,"I","D"), (2,0,1,0,0,"D","I"), (3,0,1,0,0,"D","I"),
    (4,0,2,0,0,"D","I"), (4,1,5,0,4,"I","D"), (4,1,5,1,4,"I","D"),
    (4,1,6,0,4,"D","D"), (5,0,6,0,5,"D","D"), (6,0,5,0,5,"D","D"),
    (6,0,5,1,5,"I","D"), (6,0,6,0,6,"D","I"), (6,0,6,0,6,"I","D"),
)


def verify():
    replay = Path(__file__).resolve().parent
    setup = run_path(str(replay/"ch39-minimum-orbit-specific-probability.py"))["setup"]
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    optimal = setup()
    support = [("class",v) for v in CLASSES]
    n = 6636
    offsets = (n-2345666,n-964890)
    assert optimal(offsets,support).is_zero()
    radius = arb("1e-11")
    assert optimal(offsets,support,radius).is_zero()
    kappa = arb(23)/25*arb(2)**-20
    shift = shell["soft"]["point"]["study"]["offset"]
    delta = kappa/(acb.zeta_zero(n).imag**2+shift)
    tail = kappa/(arb(10000)**2+shift)
    # 同一个零点的未乘荷钟相位同时逼近参考值，两荷的相位误差至多乘四。
    assert 0 < delta and 0 < tail
    assert 18*4*arb("1e-14")+6*(delta+tail)+arb("1e-46") < radius
    # 两个正相对宽度的高度窗也可满足指定中央六能量的必要高度比例。
    assert arb("2.26")/arb("1.001") > arb("2.25390") and arb("2.26226") < arb("2.26860")
    assert (arb("2.26")-arb("1.001"))*600000 > arb("728522.42")
    # 对具体 n 再检验全域唯一性。鲁棒邻域只断言目标多项式的实现成本。
    selector = run_path(str(replay/"ch39-rational-unique-light-action.py"))
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    reader = run_path(str(replay/"ch39-light-spectrum-zero-selector.py"))
    light = reader["light_spectra"](replay)
    ctx.prec = 320
    actual = [v for q in (-4,2) for v in light(q,n)]
    ns = [(e*10**46+arb(1)/2).floor().unique_fmpz() for e in actual]
    assert all(k is not None for k in ns)
    roots = [arb(k)/10**46 for k in ns]
    eps = arb("1e-32")
    selector["verify_number_boxes"](pump,shell,((-4,n),(2,n)),[x+arb(0,eps.upper()) for x in roots])
    lower = eps**2*min(prod((abs(x-y)/2)**2 for j,y in enumerate(roots) if i!=j).lower()
                       for i,x in enumerate(roots))
    upper = sum(prod((x-y)**2 for y in roots) for x in actual)
    assert lower-upper > arb("7.85e-69")
    print("其他实际编号相对目标的统一作用差下界", (lower-upper).str(12))
    print("两荷编号、具名轨道费用及独立参考根半宽", n, 0, radius.str(4))


if __name__ == "__main__":
    verify()
