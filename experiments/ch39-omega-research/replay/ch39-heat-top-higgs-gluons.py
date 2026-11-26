"""顶层热迹斜率约束、希格斯胶子振幅及任意轻门户作补偿的分宽下界。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


broken = run_path(str(Path(__file__).with_name("ch39-broken-electroweak-vacuum.py")))
ctx.prec = 224
mh, v = arb("125.25"), arb("246.22")


def loop(mass):
    x = (2*mass/mh)**2
    if x > 1:
        f = (1/x.sqrt()).asin()**2
    else:
        # 用 (1+β)²/x 避免 μ 层在 1−β 中消位。
        z = ((1+(1-x).sqrt())**2/x).log()
        f = -acb(z, -arb.pi())**2/4
    return 2*x*(1+(1-x)*f), x*(x*f-1)


def verify():
    names = ("u","d","s","c","b")
    masses = {k:z.sqrt() for k,_,_,z in broken["light_spectrum"](arb(1)) if k in (*names,"t")}
    sm = sum(loop(m)[0] for m in masses.values())
    top = 5*loop(masses["t"])[1]
    print("标准模型胶子振幅及色六维标量的额外振幅", sm.str(9), top.str(9))
    # 每种表示给出 2 n T(R) 和实标量分量数，n=1（复）或 1/2（实）。
    representations = {"基本":(1,6), "六维":(5,12), "伴随":(3,8)}
    options = [[(n, "基本" if n==a else "六维" if n==b else "伴随", arb(1))
                for n in names] for a,b in permutations(names,2)]
    masses["差分"] = arb("5.4600e-12").union(arb("5.4614e-12"))
    for omit in (*names,"差分"):
        # 差分层 μ(r)=r²/M，其对数质量斜率是夸克层的两倍。
        options.append([(n,"伴随", arb(2) if n=="差分" else arb(1))
                        for n in (*names,"差分") if n!=omit])
    assert len(options) == 26
    amps, widths, tolerance_widths = [], [], []
    for config in options:
        fixed_light = acb(0)
        dual_norm = arb(0)
        for name,rep,log_slope in config:
            mass = masses[name]
            index, dim = representations[rep]
            scalar = loop(mass)[1]
            fixed_light += index*log_slope*scalar
            coefficient = index*v*scalar/(2*mass**2)
            phase = (1-(2*mass/mh)**2).sqrt()
            dual_norm += abs(coefficient)**2/(dim*phase/(32*arb.pi()*mh))
        amps.append(abs(sm+top+fixed_light)/abs(sm))
        for cap, results in [(arb(1), widths), (arb("1.10"), tolerance_widths)]:
            needed = abs(sm+top)-cap*abs(sm)
            assert needed>0
            results.append(1000*needed**2/dual_norm)
    low, high = min(x.lower() for x in amps), max(x.upper() for x in amps)
    min_width = min(x.lower() for x in widths)
    min_tolerant = min(x.lower() for x in tolerance_widths)
    print("26 种跟随质量时的胶子振幅倍率", low.str(9),high.str(9))
    print("只调轻门户，使振幅至多为标准值或 1.10 倍的分宽下界/MeV",
          min_width.str(9),min_tolerant.str(9))
    assert 2<low<high<3 and min_tolerant>10


if __name__ == "__main__":
    verify()
