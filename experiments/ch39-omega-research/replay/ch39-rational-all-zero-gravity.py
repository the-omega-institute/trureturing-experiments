"""两个有理截止与整数权，对原全双钟相位同时给出引力两项的相对误差。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    spectrum = run_path(str(root / "ch39-monotone-chiral-phase-obstruction.py"))
    model, *_ = spectrum["torus"]["local"]["inputs"]()
    ctx.prec = 224
    bands = spectrum["shell_bands"](model)
    assert all(h > 0 for values in bands.values() for h in values)
    # 质量比例均为正，先取壳带下端再平方可保留全球的严格质量间隙。
    lower = {q: [h.lower() for h in values] for q, values in bands.items()}
    fields = run_path(str(root / "ch39-chiral-minimal-neutral-gravity.py"))["fields"](lower)
    floor = min(z.lower() for _, _, z in fields[1:])
    assert len(fields) == 99 and floor > 0
    ctx.prec = 640
    epsilon = arb(10)**-122
    t = epsilon/93
    s = t/arb(10)**130
    ideal = 1116*arb.pi()/epsilon
    height = (ideal+arb(1)/2).floor().unique_fmpz()
    assert height is not None and height > 1
    # (1-height) 在 s 处，height 在 t 处。两个主矩直接按非负轮廓积分。
    light_k = (s + height*(t-s))/(12*arb.pi())
    light_w = 31*(s*s + height*(t*t-s*s))/(4*arb.pi())
    # K' 和 W' 的重层相对各自无质量导数的最大偏差。
    tail = factorial(4)*(t/floor)**4
    dk = 6*sum(abs(a) for _, a, _ in fields[1:])*tail
    dw = arb(sum(abs(g) for g, _, _ in fields[1:]))*tail/62
    k = light_k*(1+arb(0,dk.upper()))
    w = light_w*(1+arb(0,dw.upper()))
    errors = [abs(k-1), abs(w/epsilon-1), abs(w/(epsilon*k*k)-1)]
    assert all(error < arb("2e-126") for error in errors)
    print("整数窗高", height)
    print("全相位平方质量下界，K、W、ΛG 的相对误差上界",
          floor.str(12), *[error.upper().str(12) for error in errors])


if __name__ == "__main__":
    verify()
