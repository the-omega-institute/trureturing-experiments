"""破缺背景上非负累计轮廓匹配真空项时，最低截止能量的锐上限。"""
from pathlib import Path
from runpy import run_path
from flint import arb, ctx


def verify(root):
    old = run_path(str(root/"ch39-chiral-induced-gravity-vacuum.py"))
    spectrum = old["spectrum"]
    background, *_ = spectrum["torus"]["local"]["inputs"]()
    ctx.prec = 224
    assert old["prefix_cap"](spectrum["shell_bands"](background),32) == 0
    index=2_306_888_100_085
    _,_,groups=run_path(str(root/"ch39-two-zero-cutoff-gravity.py"))["gravity_fields"](
        root,(index-1,index,index+1))
    ctx.prec=256
    r=(arb("246.22")/arb("1.22089e19"))**2/2
    gw=arb(".63811")
    light=[r*k for k in (gw**2/2,(gw**2+3*arb(".467453")**2/5)/2,arb(".8616")**2)]
    scale=(8*arb.pi()*arb("1e-122")/59).sqrt()
    precision=arb("1e-30")
    hi=scale*(1+precision)
    lo=scale*(1-precision)
    for number,group in zip((index-1,index,index+1),groups,strict=True):
        z0=min(z.lower() for _,_,z in group[1:])
        # 原十六双重态的每个实分量可改平方质量至多 2r。
        assert 2*r<z0
        extra=64*2*r/(arb(1).exp()*(z0-2*r))
        assert extra<arb("6e-22")
        # 全谱除零质量轻行之外，均在此公共质量下限之上。
        floor=min(z0-2*r,*light)
        assert floor>0
        # 原重谱的绝对偏振数加上分出十六个复双重态所需的保守修正。
        total=sum(abs(g) for g,_,_ in group[1:])+128+27
        tail=total*(-floor/hi).exp()
        assert tail < arb("1e-50")
        def bracket(x):
            return x*x*(59+arb(0,tail.upper()))/(8*arb.pi())
        assert bracket(lo)<arb("1e-122")<bracket(hi)
        print("零点及重移质给统一超迹上界的修正",number,extra.upper().str(8))
    mev=arb("1.22089e19")*lo.union(hi).sqrt()*arb("1e12")
    print("锐上限的截止能量（毫电子伏）",mev.str(21))
    assert arb("3.11905834380287") < mev < arb("3.11905834380289")


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
