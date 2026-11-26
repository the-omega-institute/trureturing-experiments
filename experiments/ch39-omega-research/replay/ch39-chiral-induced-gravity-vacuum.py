"""三代规范谱在全双钟相位的诱导牛顿上界和无量纲宇宙学常数下界。"""
from collections import defaultdict
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).parent
spectrum = run_path(str(ROOT / "ch39-monotone-chiral-phase-obstruction.py"))


def prefix_cap(bands, ratio):
    # 质量单位取 Rκ。偶模归一质量固定，奇模随 R 增大只会降低。
    events = defaultdict(int)
    counts = {1:-24,4:-12,2:-12,3:-4,6:-4}
    for charge,signed in counts.items():
        for h in bands[charge][3:]:
            upper = h.upper()/(ratio if charge%2 else 1)
            events[int((100*upper).ceil().unique_fmpz())] += signed
    for h in bands[0]:
        events[int((100*h.lower()).floor().unique_fmpz())] += 36
    current = maximum = 0
    for mass in sorted(events):
        current += events[mass]
        maximum = max(maximum,current)
    assert current==-248
    return maximum


def newton_fraction(neutral, ceiling):
    # 用中性质量平方的上端替代实际值，为负矢量项给出统一上包络。
    a = [v.upper()**2 for v in neutral]
    def slope(x):
        return ceiling-6*sum((-m/x).exp() for m in a)
    lo,hi = arb(1),arb(1000000)
    assert slope(lo)>0>slope(hi)
    for _ in range(48):
        mid = (lo+hi)/2
        if slope(mid)>0:
            lo=mid
        else:
            assert slope(mid)<0
            hi=mid
    x = lo.union(hi)
    least = sum((m/x.lower()).expint(2) for m in a)
    bound = (x*(ceiling-6*least)/(2*arb.pi())).upper()
    unit = arb("2.435e18")/724**2
    fraction = bound*(32*unit/arb("1.22089e19"))**2
    return fraction


def verify():
    model,_,_,_,_ = spectrum["torus"]["local"]["inputs"]()
    ctx.prec = 192
    bands = spectrum["shell_bands"](model)
    for ratio,target in ((4,16),(16,0)):
        cap = prefix_cap(bands,ratio)
        assert cap==target
        for full in (False,True):
            # 至多另外三个复双重态，各四个实偏振。
            deficit = 62-cap-(12 if full else 0)
            # 牛顿曲率权含轻层 45 Weyl、一个复双重态、十二轻规范矢量；
            # 重层 240 Dirac 分量、十六或十九个复双重态与十八层十二个 Proca。
            ceiling = arb(45)/6+arb(4)/6-8+arb(240)/3+arb(4*(19 if full else 16))/6
            dimless = deficit*arb.pi()/(2*ceiling**2)
            print("倍率至少",ratio,"十九重标量" if full else "十六重标量",
                  "全谱前缀偏振亏量下界与 ΛG 下界",deficit,dimless.lower().str(11))
            assert dimless>arb(".005")
            if ratio==16:
                fraction = newton_fraction(bands[0],ceiling)
                print("固定中心倍率 32 时正牛顿倒数的观测份额上界",fraction.upper().str(10))
                assert fraction<arb(".0001")


if __name__=="__main__":
    verify()
