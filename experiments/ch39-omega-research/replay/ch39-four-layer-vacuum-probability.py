"""四层真空标量行列式作为零点赋权的有限三点检验。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpz_poly

physical = run_path(str(Path(__file__).with_name("ch39-first-three-mode-observables.py")))


def verify():
    ctx.prec, ctx.cap = 320, 3
    vacuum = physical["joint"]["vacuum"]
    count = 4
    low_rows = vacuum["low_sector"](count, arb(1))
    assert [z for _, z in low_rows] == [0, 1, 2, 3]
    # 消去四层以外已逐层配平的净偏振。四层净偏振除以 4 为负三阶差分系数。
    t = fmpz_poly([0, 1])
    numerator = (1+t)**3*(1+3*t)
    denominator = (1+2*t)**3
    assert numerator-denominator == 2*t**3+3*t**4
    assert t*(numerator.derivative()*denominator-numerator*denominator.derivative())*(1+t)*(1+2*t)*(1+3*t) == 6*t**3*numerator*denominator

    coefficient = (9*arb(3).log()-12*arb(2).log())/2
    mass_floor = 32*arb("2.435e18")*9/724**2
    weak_cap = (arb("246.22")**2/(mass_floor*arb("1.22089e19")))**2
    error = weak_cap/(arb(".03")*coefficient)
    assert 0 < error < arb("1e-54")
    lower = (arb(".998")*arb.pi()*arb("1e-122")/coefficient).sqrt()
    upper = (arb("1.002")*arb.pi()*arb("1e-122")/(coefficient*(1-error))).sqrt()
    spread = upper/lower
    assert 1 < spread < arb("1.00201")
    first, second = [acb.zeta_zero(n).imag for n in (8359, 12336)]
    bound = 2*(second/first).log()+spread.log()
    exponent = arb(259).log()/bound
    assert exponent > 8
    four_layer = (2*bound).exp()
    assert four_layer < 4
    t_max = upper/arb(".00015")
    at_shared_scale = arb(259).log()/(2*t_max**3+3*t_max**4)
    assert at_shared_scale > arb("1.16746e171")
    print("四层质量平方比上界、单个高斯四层模型的三点校正上界", spread.str(13), four_layer.str(13))
    print("自由核尺度所需绝对幂下界", exponent.str(13))
    print("共同核尺度所需绝对幂下界", at_shared_scale.str(13))


if __name__ == "__main__":
    verify()
