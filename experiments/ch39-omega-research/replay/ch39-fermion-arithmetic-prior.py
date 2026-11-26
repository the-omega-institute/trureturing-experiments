"""零点参考态先验与单态乘法链的必要刚度和充分温度窗。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
angle = run_path(str(root/"ch39-first-arithmetic-plaquette-winding.py"))["angle"]
addition_chain = run_path(str(root/"ch39-neutral-fermion-phase-register.py"))["addition_chain"]


def energy(offset):
    # 取 M=1,r=1/2。相位无关的加法常数在所有比较中消去。
    h = [arb(5)/4+((offset+2*k*arb.pi())/5).cos() for k in range(5)]
    return -2*sum(x**2*x.log() for x in h)/(16*arb.pi()**2)


def verify():
    ctx.prec = 180
    power = 23226
    _, weights = addition_chain(power)
    weights = weights[1:]
    assert sum(weights) == power-1
    eps_target = abs(power*angle(11431)-2*arb.pi()/3)
    eps_second = abs(power*angle(1144)-14*arb.pi()/3)
    inverse_weight = lambda n: arb(1)/4+acb.zeta_zero(n).imag**2
    odds = inverse_weight(11431)/inverse_weight(1144)
    curvature = 2*odds.log()/(eps_second**2-eps_target**2)
    equal_needed = curvature*sum(x*x for x in weights)
    any_needed = curvature*(power-1)**2
    assert equal_needed > arb("2.8309e17") and any_needed > arb("1.3790e18")

    kappa = arb("1e23")
    inverse_temperature = arb("1e15")
    margin = arb("2e-6")
    single = margin/sum(weights)
    # 直接用半角平方避免在大刚度下相消。
    tail = (len(weights)*arb.pi()*kappa.sqrt()*
            (arb(1)/2-2*kappa*(single/2).sin()**2).exp())
    gap = energy(eps_second-margin)-energy(eps_target+margin)
    bias = energy(eps_target+margin)-energy(arb(0))
    total = 1+arb.const_euler()/2-(4*arb.pi()).log()/2
    error = (2*total*inverse_weight(11431)*
             ((-inverse_temperature*gap).exp()+tail*(inverse_temperature*bias).exp())/(1-tail))
    assert gap > arb("1.3532e-13") and tail < arb("1e-145")
    assert error < arb("1e-51")
    print("共同刚度必要下界",equal_needed.str(12),"任意正乘法链总刚度下界",
          any_needed.str(12),"有限前缀错配概率上界",error.str(12))


if __name__ == "__main__":
    verify()
