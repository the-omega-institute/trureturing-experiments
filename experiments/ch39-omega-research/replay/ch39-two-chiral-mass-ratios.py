"""用右弱壳两个稳定取向消去未知端口强度，从三项矢量质量恢复零点高度。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

original = run_path(str(Path(__file__).with_name("ch39-zero-biased-weak-plaquette.py")))


def coefficients(height):
    den = height**2+arb(1)/4
    a,b = (height**2-arb(1)/4)/den,height/den
    root = arb(3).sqrt()
    return 2+a-root*b,2+a+root*b,b/(2+a)


def verify():
    original["geometry"]()
    ctx.prec = 160
    assert acb.zeta_zero(1).imag>14
    root = arb(3).sqrt()
    # 归一 g=1，kappa/g>=2，lambda/g>=4。适用于全部正高度。
    max_energy = 3+root/14
    split = 2*root/14
    gap = 6-max_energy
    assert max_energy<8 and gap>split
    assert 4*(1-split/gap)>split
    heights = [acb.zeta_zero(n).imag for n in (11430,11431,11432)]
    data = [coefficients(h) for h in heights]
    s0 = data[1][2]
    eps = arb("1e-10")
    eta = ((1+eps)/(1-eps))**2-1
    # 三个观测质量各有独立误差，区间按更宽的两份叶质量误差处理。
    sum_min = (data[1][0]+data[1][1])/12
    factor = 1+2/sum_min
    radius = (1/root+s0)*factor*eta/(1-factor*eta)
    inferred = s0+arb(0,radius.upper())
    assert data[2][2] < inferred < data[0][2]
    recovered = (1+(1-3*inferred**2).sqrt())/(6*inferred)
    assert heights[0] < recovered < heights[2] and recovered.contains(heights[1])
    print("calibration-free shape",s0.str(14),"error radius",radius.str(10),
          "height endpoints",recovered.lower().str(12),recovered.upper().str(12),flush=True)
    for h,(plus,minus,s) in zip(heights,data,strict=True):
        assert s.overlaps(h/(3*h**2+arb(1)/4))
        for eigen,nominal in zip((plus,minus),data[1][:2],strict=True):
            retuned = arb(1)/5*nominal/eigen
            assert arb(1)/6 < retuned < arb(1)/4
    # 未重新标定的一个比值会随端口耦合改变。
    assert data[0][0] != data[1][0] and data[2][1] != data[1][1]


if __name__=="__main__":
    verify()
