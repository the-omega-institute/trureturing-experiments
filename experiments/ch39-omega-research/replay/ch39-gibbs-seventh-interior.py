"""全支撑吉布斯解保持六项时间时，第七项距两个全谱端点的双侧裕量。"""
from math import prod
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
scale=2_903_600_000
indices=(8359,349501044,142269938072,142269606336,295872524995,295668865161,13219432687559,143921293989490)


def verify():
    ctx.prec=768
    pinned=run_path(str(root/"ch39-pinned-cp-zero-support.py"))
    extra=run_path(str(root/"ch39-seventh-time-zero-dictionary.py"))
    assert set(indices) <= set(pinned["indices"]) | set(extra["spectrum_upper"])
    ctx.prec=144
    heights=sorted((acb.zeta_zero(n).imag for n in indices), key=float)
    ctx.prec=768
    y=[1/(1+(h/scale)**2) for h in heights]
    assert all(a>b for a,b in zip(y,y[1:]))
    def D(t):
        return prod(1+(j*j-1)*t for j in range(2,7))
    v=[D(t)/prod(t-s for j,s in enumerate(y) if j!=i) for i,t in enumerate(y)]
    absolute=sum(abs(w) for w in v)
    w=[t/absolute for t in v]
    assert sum(w).contains(0)
    for j in range(1,7):
        assert sum(a*(1-t)/(1+(j*j-1)*t) for t,a in zip(y,w)).contains(0)
    # 七阶差商：将 (1-y)D(y)/(1+48y) 作带余除法，至多五次多项式被消去。
    pole=-arb(1)/48
    seventh=abs((1-pole)*D(pole)*(-48)**7/(prod(1+48*t for t in y)*absolute))
    direct=abs(sum(a*(1-t)/(1+48*t) for t,a in zip(y,w)))
    assert direct.overlaps(seventh) and seventh>0
    h=arb(8000)
    prior_cap=((h/(2*arb.pi())).log()+1)/(2*arb.pi()*h)+arb(".28")*(2*h.log()+arb(".5"))/h**2
    min_prior=1/(max(z.upper() for z in heights)**2*prior_cap)
    factor=min_prior*seventh
    print("单位一范数方向的第七读数、最小先验下界、双侧系数",seventh.str(15),min_prior.str(15),factor.str(15))
    assert factor > arb("1.14198e-34")


if __name__=="__main__":
    verify()
