"""第 8359 枚在连续双钟壳上恢复参考味和三个物理窗口的最小转角。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


root = Path(__file__).parent
flavor = run_path(str(root/"ch39-flavor-invariants.py"))
physical = run_path(str(root/"ch39-first-three-mode-observables.py"))
shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
window = physical["mask"]["window"]


def verify():
    model,initial=physical["model_inputs"]()
    ctx.prec=224
    height=acb.zeta_zero(8359).imag
    original=shell["geo"]["phases"](height)
    other,delta=flavor["reference_rotations"]()
    assert abs(other)/2>arb(".12") and 0<delta<arb(".051")
    def phases(t):
        return {724:original[724]*acb(0,t+delta/2).exp(),
                804:original[804]*acb(0,t-delta/2).exp()}

    neutral=model["shell"]["shell_spectra"]((0,),zero_index=8359)[0]
    center=sum(neutral)/18
    ordered=physical["joint"]["chiral"]["ordered"]
    def read(t):
        levels={0:neutral}
        for q in model["charges"]:
            eigen=shell["laplacian"](q,height,phases=phases(t)).eig(multiple=True)
            assert len(eigen)==19 and all(z.imag.contains(0) for z in eigen)
            levels[q]=[z.real for z in eigen]
        layers={q:ordered([v*abs(v-center) for v in values]) for q,values in levels.items()}
        return physical["layer_observables"](model,initial,layers)

    scale=arb("246.22")**2/(32*arb("2.435e18")/724**2*arb("1.22089e19"))
    coefficient=(9*arb(3).log()-12*arb(2).log())/2
    boundary=scale*(coefficient/(arb.pi()*arb("1e-122")*arb("1.00200001"))).sqrt().sqrt()
    assert boundary>center**2/4
    assert arb("1.00200001")*(1-(scale/9)**2/(arb(".03")*coefficient))>arb("1.002")
    separator=(center+(center**2+4*boundary).sqrt())/2
    count=window["supplier"](2,separator)
    lower,upper=arb(".00264225"),arb(".00264228")
    # 固定单位辐边的正谱序给出实际壳的上界。区间宽导致消元枢轴未定时平分。
    stack=[(-arb(1),arb(1))]
    while stack:
        left,right=stack.pop()
        try:
            excluded=1+count(phases(lower*left.union(right)))>=16
        except AssertionError:
            excluded=False
        if not excluded:
            assert right-left>arb("1e-12")
            mid=(left+right)/2
            stack.extend(((left,mid),(mid,right)))
    for shift,kind in ((arb(0),"最近味相位"),(-upper,"进入曲率窗"),(arb("-.004840689"),"曲率近中心")):
        energy,strong,curve=read(shift)
        assert arb(".117")<strong<arb(".119")
        if shift==0:
            assert arb("1.00440")<curve<arb("1.00442")
        else:
            assert abs(curve-1)<(arb(".002") if kind=="进入曲率窗" else arb("1e-9"))
        print(kind,energy.str(8),strong.str(9),curve.str(12))
    assert arb(".02768915")<delta/2+lower<delta/2+upper<arb(".02768920")
    print("同时达标的最大单钟转角下确界在",(delta/2+lower).str(13),(delta/2+upper).str(13),"之间")


if __name__=="__main__":
    verify()
