"""接近最小共同尺度时实际零点编号及高处算术份额的下界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

floor = run_path(str(Path(__file__).with_name("ch39-neutral-five-scale-floor.py")))


def verify():
    ctx.prec, ctx.cap = 512, 3
    first = acb.zeta_zero(8359).imag
    scale_floor = first*arb(117387924448).sqrt()
    scale_cap = 2_903_600_000
    right = ((scale_cap/first)**2).upper()
    left = arb("1.455967e-8")
    margin, *_ = floor["envelope"](right, left)
    assert margin < 0
    cutoff = 24027719528722
    assert scale_floor/left.sqrt() > cutoff
    count = arb(cutoff).zeta_nzeros().unique_fmpz()
    assert count == 106969920543670
    print("必须超过的零点高度及该高度以下正零点数", cutoff, count)
    guarantees = ((5, ".37426"), (10, ".32366"), (15, ".23893"), (20, ".12026"), (23, ".03277"))
    for trillions, guarantee in guarantees:
        height = trillions*10**12
        lower = ((scale_floor/height)**2).lower()
        # x<lower 即推出 γ>height。小 x 时 ηγg 至多为 η₀g。
        margin, _, _, contact, _, _, quotient, den = floor["envelope"](right, lower)
        interval = arb(0).union(lower)
        debt = (-contact(interval)*quotient(interval)/den(interval)).upper()
        assert debt > 0
        probability = -margin/debt
        assert arb(guarantee) < probability < 1
        if trillions == 15:
            h = arb(height)
            tail = ((h/(2*arb.pi())).log()+1)/(2*arb.pi()*h)
            tail += arb(".28")*(2*h.log()+arb(".5"))/h**2
            gain = probability/(arb("4.66881e-11")*first**2*tail)
            assert gain > arb("2.28327e14")
            assert arb("1.680")**63 < gain and arb("2.147")**43 < gain
            print("相对首枚的逆平方谱效率倍率下界", gain.str(15))
        print("高度以上的算术概率下界", height, probability.str(14))


if __name__ == "__main__":
    verify()
