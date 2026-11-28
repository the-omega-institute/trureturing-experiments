"""用荷拉普拉斯的四步闭路迹，在留出重谱上预报真空矩变化方向。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).resolve().parent
SOURCE = 2_345_698
HOLDOUT = tuple(range(2_345_724, 2_345_771, 2))


def full_trace(geometry, number):
    h = acb.zeta_zero(number).imag
    offset = geometry["soft"]["point"]["study"]["offset"]
    delta = arb(23)/25*arb(2)**-20/(h*h+offset)
    pivot = 15-delta/3
    unit = arb("2.435e18")/(724**2*arb("1.22089e19"))
    total = arb(0)
    for q, g in ((0,36),(1,-24),(4,-12),(2,-12),(3,-4),(6,-4)):
        laplace = geometry["laplacian"](q,h)
        polynomial = laplace*laplace - laplace*pivot
        trace = sum(abs(v)**2 for v in polynomial.entries())
        scale = unit if q%2 else 32*unit
        total += g*scale**2*trace
    return total


def verify():
    load = run_path(str(ROOT/"ch39-four-background-crossing.py"))["load"]
    geometry = load("three-ray-chiral-boundary")
    ctx.prec = 320
    numbers = (SOURCE, *HOLDOUT)
    # 牛顿项的系数与重谱输入无关。按原场表的重数直接加和。
    alpha = arb(49)/6-8 + 16*(arb(6+3+3+2+1)/3 + arb(4)/6) - 18*6
    assert (alpha+arb(103)/6).contains(0)
    h0 = full_trace(geometry,SOURCE)
    prediction = {n:(full_trace(geometry,n)-h0)/(2*alpha) for n in HOLDOUT}
    hits, wrong, open_cases = [], [], []
    rising, falling = [], []
    order_check = load("fifth-shell-vacuum-order")
    for a,b,lead,change,actual in order_check["vacuum_intervals"](load, numbers, source_only=True):
        guess = prediction[b]
        assert not guess.contains(0)
        (rising if guess>0 else falling).append(b)
        sign = guess*change
        (hits if sign>0 else wrong if sign<0 else open_cases).append(b)
        print(b, "全谱预报、重谱预报、真空矩差下上端、冻结见证", guess.str(9), lead.str(9),
              change.lower().str(9), change.upper().str(9), actual.str(9), flush=True)
    assert hits == list(HOLDOUT) and not wrong and not open_cases
    print("真空矩升高、降低的编号",rising,falling,flush=True)


if __name__=="__main__":
    verify()
