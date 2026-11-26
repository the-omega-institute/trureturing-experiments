"""Cayley 质量桥任意尺度和左右斜率分支上的胶子与顶伴随产生比较。"""
from itertools import permutations, product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


higgs = run_path(str(Path(__file__).with_name("ch39-heat-top-higgs-gluons.py")))


def verify():
    ctx.prec = 224
    names = ("u","d","s","c","b")
    masses = {q:z.sqrt() for q,_,_,z in higgs["broken"]["light_spectrum"](arb(1)) if q in (*names,"t")}
    forms = {q:higgs["loop"](m) for q,m in masses.items()}
    sm = sum(f[0] for f in forms.values())
    top = sum(x*w for x,w in zip(forms["t"],(1,5)))
    assert top > 0
    slopes = {q:(1-(masses[q]/masses["t"])**2).sqrt().lower().union(arb(1)) for q in names}
    ultra = 6*higgs["loop"](arb("5.4600e-12").union(arb("5.4614e-12")))[1]
    options = [({q:(1 if q==a else 5 if q==b else 3) for q in names}, False)
               for a,b in permutations(names,2)]
    options.extend(({q:(0 if q==omit else 3) for q in names}, omit!="差分") for omit in (*names,"差分"))
    bounds, near_top = [], []
    for indices, use_ultra in options:
        blocks = [slopes[q]*(forms[q][0]+indices[q]*forms[q][1]) for q in names]
        for signs in product((-1,1), repeat=5):
            rem = sum((s*b for s,b in zip(signs,blocks)), acb(ultra) if use_ultra else acb(0))
            assert abs(rem.imag)>0
            # 任意非零实 κ_t，|A+B/κ_t| ≥ 实 A 到实直线 span(B) 的距离。
            bounds.append((top*abs(rem.imag)/(abs(rem)*abs(sm)))**2)
            # 若 |κ_t|≥0.8，反三角不等式直接给出较强界。
            hard = top-abs(rem)/arb(".8")
            assert hard>0
            near_top.append((hard/abs(sm))**2)
    print("两种有理斜率分支、26 种分配，通用及 |κ_t|≥0.8 时归一产生率比下界",
          min(b.lower() for b in bounds).str(9), min(b.lower() for b in near_top).str(9))
    assert len(bounds)==832
    assert min(b.lower() for b in bounds)>1


if __name__=="__main__":
    verify()
