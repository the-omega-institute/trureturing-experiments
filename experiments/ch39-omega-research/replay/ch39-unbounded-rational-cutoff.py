"""每个请求位数只用四个严格正的有理截止，全钟域逼近两个引力目标。"""
from argparse import ArgumentParser
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def dyadic_pi(m):
    # Machin 恒等式，两项反正切均用有限交错和。
    def partial(q):
        return sum(Fraction((-1)**j, (2*j+1)*q**(2*j+1)) for j in range(m))

    scale = 1 << m
    value = 16*partial(5)-4*partial(239)
    remainder = Fraction(20, (2*m+1)*5**(2*m+1))
    assert remainder*scale < 1
    return Fraction((value*scale).__floor__(), scale)


def verify(bits):
    root = Path(__file__).resolve().parent
    spectrum = run_path(str(root / "ch39-monotone-chiral-phase-obstruction.py"))
    model, *_ = spectrum["torus"]["local"]["inputs"]()
    ctx.prec = 224
    bands = spectrum["shell_bands"](model)
    assert all(h > 0 for values in bands.values() for h in values)
    fields = run_path(str(root / "ch39-chiral-minimal-neutral-gravity.py"))["fields"]
    low = fields({q: [h.lower() for h in values] for q, values in bands.items()})
    high = fields({q: [h.upper() for h in values] for q, values in bands.items()})
    assert min(z.lower() for _, _, z in low[1:]) > arb("4.55e-12")
    assert max(z.upper() for _, _, z in high) < 1
    alpha = sum(abs(a) for _, a, _ in high)
    spin = sum(abs(g) for g, _, _ in high)
    assert alpha < 2000 and spin < 2000 and sum(g for g, _, _ in high) == -310

    # 以下界均为相对误差除以 s=2^-m 的上界，对一切 m≥512 有效。
    pi = arb.pi()
    epsilon = arb(10)**-122
    k_ir = 2/pi + 1/(12*pi) + 49*alpha/(pi*arb("4.55e-12"))
    k_uv = 16*epsilon*alpha/(465*2*pi)
    w_ir = 145*spin/(8*pi*epsilon)
    w_uv = 2/pi + 4*spin/(310*pi)
    cap = max(k_ir+k_uv, w_ir+w_uv)
    assert cap < arb(2)**420
    # 若 |K-1|、|W/ε-1|≤e≤1/4，曲率的相对误差至多 6e。
    assert (3+arb(1)/4)/(1-arb(1)/4)**2 < 6

    m = bits + 512
    T = 1 << m
    s = Fraction(1, T)
    p = dyadic_pi(m)
    assert 3 < p < 4
    A, B = 12*p*T, Fraction(4, 465*10**122*T*T)*p
    assert 1 < A < 48*T and 0 < B < 1
    assert s < 2*s < T < 2*T
    print("所求相对精度位数、截止指数 m、统一误差系数的上界",
          bits, m, cap.upper().str(12))


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--bits", type=int, default=256)
    args = parser.parse_args()
    if args.bits < 0:
        parser.error("精度位数须为非负整数")
    verify(args.bits)
