"""参考混合的行列式尺度与五荷三代共用的算术前缀。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


source=run_path(str(Path(__file__).with_name("ch39-flavor-scale-time-prefix.py")))
mask=source["clock"]["mask"]
indices=(2556787,2556892,3936059,3936275,4767305,4767415,6197447,6197558,
         8758142,8758272,15269054,15269427,55329118,55329394)


def verify():
    ctx.prec,ctx.cap=1536,3
    times=source["clock"]["clock"]["times"]()
    lo,hi=source["determinant"]["scale_interval"]()
    scale=lo.union(hi)
    ctx.prec=256
    heights=[acb.zeta_zero(n).imag for n in indices]
    accepted=mask["three_mode_test"]()
    assert all(accepted(h) for h in heights)
    first,last=indices[:2]
    gap=[z.imag for z in acb.zeta_zeros(first+1,last-first-1)]
    assert len(gap)==last-first-1 and all(not accepted(h) for h in gap)
    assert heights[0]<scale*acb.zeta_zero(1).imag<heights[1]
    ctx.prec=1536
    weights=source["certify_weights"](times,scale,13,heights)
    margin=source["separation_margin"](times,scale,14,0,heights[0],heights[1])
    assert margin>arb("6.10e-25")
    print("五荷三模十四枚的最小原概率、十四项最大误差下界",min(v.lower() for v in weights).str(12),margin.str(12))


if __name__=="__main__":
    verify()
