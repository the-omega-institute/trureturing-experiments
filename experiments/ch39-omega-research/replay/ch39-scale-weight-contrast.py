"""四十一项时间读数要求的大尺度谱权上升与温和调制的排除。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_poly, ctx


local = run_path(str(Path(__file__).with_name("ch39-three-mode-time-localization.py")))


def verify():
    ctx.prec,ctx.cap = 1536,3
    times = local["clock"]["times"]()
    moments = local["clock"]["cleared_moments"](times,41)
    poles = local["partial_fractions"](41)
    height = acb.zeta_zero(1).imag
    x = 200/height**2
    p = local["window_polynomial"](moments,x)
    eps = arb("4e-16")
    radius = eps*x
    q = arb_poly([1])
    for rate,_ in poles:
        q *= arb_poly([1,arb(rate)/200])
    numerator = p**2*(1-(1/radius**2)*arb_poly([-x,1])**2)
    residues = local["residues_of"](numerator,poles)
    signal = sum(c*t for c,t in zip(residues,times[1:],strict=True))
    signal -= arb(2)**-176*sum(abs(c) for c in residues)
    offset = arb(0,radius)
    maximum = (p(arb_poly([x,1]))(offset)**2/q(x+offset)).upper()
    mass = signal/maximum
    assert mass > arb(".16895")

    # 开窗口严格落在 (T,(1+s)T) 内，从而两端的零点计数约定均可使用。
    left = 1-arb("2.001e-16")
    s = arb("4.003e-16")
    assert left < 1/(1+eps).sqrt() < 1/(1-eps).sqrt() < left*(1+s)
    assert left*height*arb("588.4188") > 8300
    split = arb(10)**22
    two_pi = 2*arb.pi()
    threshold = arb("4.7e14")
    # 界 |N(t)-L(t)| <= 0.28 log(t) 引自 Brent--Platt--Trudgian。
    upper = split*(1+s)
    count_cap = upper/two_pi*((upper/two_pi).log()-1)+arb(7)/8+arb(".28")*upper.log()
    assert count_cap < arb(10)**23
    assert mass/count_cap > threshold*arb("7.28e-51")

    assert split.log()*arb(".06") > two_pi.log()+1
    assert (1+s).log() < arb(".001")*split.log()
    # 用 log(T) 约去上下界的共同因子，再统一用 T >= split 包住 1/T。
    ratio = (arb("1.001")*s/two_pi+arb(".561")/split)/(arb(".94")/two_pi-arb(".28")/split)
    assert ratio>0 and mass/(1-mass) > threshold*ratio
    print("第一窄窗最低质量",mass.str(13),"高处窄窗与此前计数的比值上界",ratio.str(12))
    print("所有大尺度均须有后点权重超过前点权重的",threshold.str(4),"倍")


if __name__ == "__main__":
    verify()
