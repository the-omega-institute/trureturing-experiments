"""不指定谱目标时三支规范流的最低输入改动，以及同一组输入的两份谱时间表。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
flow = run_path(str(root/"ch39-floating-scale-winner.py"))
ratios = flow["ratios"]
tolerance = arb(2)**-21


def gap_at_split(ppm, split):
    radius = arb(ppm)/10**6
    inverse = [(1/(g*(1+radius))**2).union(1/(g*(1-radius))**2) for g in flow["centers"]]
    logs = [2*(y*arb(".99").union(arb("1.01"))).log() for y in flow["yukawa_centers"]]
    hyper, weak, color = ratios
    projected = inverse[0]/hyper - inverse[1]/weak
    running, loop = flow["running"], flow["loop"]
    linear = running["BETA"][0]/hyper - running["BETA"][1]/weak
    def contrast(key):
        return [a/hyper - b/weak for a,b in zip(*running[key][:2])]
    gauge_row, yukawa_row = contrast("B"), contrast("D")

    def accumulate(step, gauge, traces):
        nonlocal projected
        correction = sum(b*g - d*y for b,d,g,y in zip(gauge_row, yukawa_row, gauge, traces))
        projected -= step*2*(linear + correction/loop)/loop

    _, _, weak_color = flow["evolve"](inverse, logs, split, 8192, observe=accumulate)
    # 任一允许的原始阈值都受 H 约束，所有重质量至少是 30 exp(27) TeV。
    mass_floor = arb(30)
    z = ((split-27).exp()/mass_floor)**2
    upper = (1+z/4).log()/(4*arb.pi()**2)
    return (projected-upper/hyper).lower(), (-weak_color-upper/color).lower()


def exclude_all_targets():
    # 在整个输入盒、整个尺度窗上重算裸规范支的导数界。
    running, loop = flow["running"], flow["loop"]
    wide = arb(".999").union(arb("1.001"))
    inverse = [1/(g*wide)**2 for g in flow["centers"]]
    logs = [2*(y*arb(".99").union(arb("1.01"))).log() for y in flow["yukawa_centers"]]
    steps = 0
    def slopes(step, gauge, traces):
        nonlocal steps
        steps += 1
        if steps <= 23*128:
            return
        rates = [-2*(beta+sum(b*g-d*y for b,d,g,y in zip(row,src,gauge,traces))/loop)/(loop*v)
                 for beta,row,src,v in zip(running["BETA"],running["B"],running["D"],ratios)]
        assert rates[0] < rates[1]
        assert rates[2]-rates[1] > 1/(2*arb.pi()**2*ratios[2])
    flow["evolve"](inverse,logs,arb(33),33*128,observe=slopes)
    assert steps == 33*128
    for ppm, text, lower in ((0,"26.8479",".0057"),(620,"26.81523",".0000169")):
        before, after = gap_at_split(ppm, arb(text))
        assert before > arb(lower) and after > arb(lower) > 2*tolerance
        print("Input ppm",ppm,"projected separation lower bounds",before.str(13),after.str(13))


def verify_readouts():
    grid = run_path(str(root/"ch39-grid-gauge-readout.py"))
    mixing = run_path(str(root/"ch39-ckm-gauge-displacement.py"))
    control = run_path(str(root/"ch39-center-free-three-pulses.py"))
    # 三个输入和匹配尺度是拟合得到的有限十进制数。九项 Yukawa 移到各自下沿。
    ppm = list(map(arb, ("594.636366", "-647.332520", "647.332520")))
    inputs = [(g*(1+s/10**6))**2 for g,s in zip(flow["centers"], ppm)]
    inputs += [(arb(".99")*y)**2 for y in flow["yukawa_centers"]]
    assert all(abs(s) < arb("647.333") for s in ppm)
    scale = arb("26.81678585294626")
    for _ in range(128):
        inputs = grid["advance"](inputs, scale/128, 14)
    drift = [a+b for a,b in mixing["displacement"](27)]
    reference = arb("26.81676714125598")
    mass_ratio = 42*(reference-27).exp()
    mass = mass_ratio*arb(".99992").union(arb("1.00008"))
    assert mass > 30

    def kernel(momentum):
        z = (momentum/mass)**2
        return (z/5-3*z**2/140).lower().union((z/5).upper())
    shift = (kernel((scale-27).exp())-kernel(arb(-27).exp()))/(4*arb.pi()**2)
    absolute = [1/(v*u)+arb(0,d.upper())-shift/v for v,u,d in zip(ratios,inputs,drift)]
    examples = (
        ((11431,), (1,)),
        ((11415,11430,11432,11389,11448),
         (20547934198306,118869088944549,850933240432901,4716893508610,4932842915634)),
    )
    ordered, study, box = (control[k] for k in ("ordered", "study", "box"))
    for labels, counts in examples:
        weights = [arb(n)/sum(counts) for n in counts]
        target = sum(w*flow["target"](label,scale) for w,label in zip(weights,labels))
        errors = [abs(v-target) for v in absolute]
        assert all(e < tolerance for e in errors)
        blocks = [study["native"]["plane_response"](acb.zeta_zero(n).imag) for n in labels]
        largest = arb(0)
        half = arb(2)**-29
        for signs in product((-1,1),repeat=4):
            cube = [v+s*half+arb(0,half) for v,s in zip(box["centers"],signs)]
            steps = {n:ordered["evolution"](block,w,cube) for n,block,w in zip(labels,blocks,weights)}
            obs = ordered["readings"](steps,labels[::-1])
            for value,center,width in zip(obs,ordered["flavor"]["shell"]["centers"],[arb(".000055")]*3+[arb(".006")]):
                largest = largest.max(abs(value-center)/width)
        assert largest < 1
        print("Labels",labels,"gauge residual uppers",*(e.upper().str(10) for e in errors),
              "maximum normalized CKM error",largest.upper().str(10))


if __name__ == "__main__":
    ctx.prec = 160
    exclude_all_targets()
    verify_readouts()
