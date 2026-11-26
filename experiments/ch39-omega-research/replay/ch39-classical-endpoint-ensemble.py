"""六条右弱环边的完整经典热积分、取向浓度与局部翻转瓶颈。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent
q = arb
heat_min, heat_max = q("1.25e-18"), q("2e-17")


def volume_bounds():
    # 归一平面测度 d²z/π。把温度上界加倍以估计受限积分。
    a_min = q(8)/11
    thermal_cap = 2*heat_max
    lo, hi, sine = q(".98"), q("1.02"), q(".7")
    cos_lo = ((lo-q(".1"))**2-2*hi**2)/(2*lo**2)
    cos_hi = ((hi+q(".1"))**2-2*lo**2)/(2*hi**2)
    assert cos_hi < 0 and 1-cos_lo**2 > sine**2
    # 半径规则区的两支相位坐标，及 |sum P| 超出小圆盘的余项。
    cutoff = q(".1")*a_min
    jac_min = (lo*a_min)**2*sine
    # exp(-u) <= 1/u²，故余项除以 H/k 后可以统一压小。
    phase_tail = 4*(thermal_cap/2)/cutoff**4
    assert phase_tail < q(".001")
    regular = (q.pi()/q("2.5"))**3*(1/(q.pi()*jac_min)+phase_tail)
    # 径向非常规区分出一半高斯衰减，余下部分对最大一条腿做复高斯积分。
    bad_radial_cost = q("1.25")*(q(".02")*a_min)**2
    crude = 2*4*(8+q(3)/2)
    exceptional = crude*(-bad_radial_cost/thermal_cap).exp()/heat_min**3
    assert exceptional < q(".001")
    upper = regular+exceptional

    # 一支取向上同时改变六个模长平方和闭合和的两个实分量。
    width, disk = q(".09"), q(".15")
    assert width*thermal_cap.sqrt() < q(".02")*a_min
    assert disk*(thermal_cap/2).sqrt() < cutoff
    path_error = width+(3*width+disk/q(2).sqrt())/sine
    penalty = 45*width**2+q(5)/3*disk**2+4*(width**2+2*path_error**2)
    lower = (2*width)**6*disk**2/(4*q.pi()*hi**2)*(-penalty).exp()
    spread = (upper/lower).log()
    assert spread < 21
    print("去除最低能后，配分因子除以 H^4/k 的上下界",lower.str(10),upper.str(10),
          "对数跨度上界",spread.str(10),flush=True)
    return spread


def verify():
    ctx.prec = 128
    spread = volume_bounds()
    # 既有平衡证书给非目标总权概率界，重加权最多乘 exp(spread)。
    balance = run_path(str(root/"ch39-thermal-plaquette-occupation.py"))
    balance["verify"]()
    assert 4000-spread > 3979 and 199984-spread > 199963

    # 动力学只需更改与端点扇区有关的能差；尺度纯步仍有原来的 190 H。
    budgets = [(q(200000)/8,10001),(q(800000)/8,50001),(q(1000000)/8,24385),
               (q(299000000)*q("2e-8")/(8*q(".0005")),661)]
    count_risk = sum((-(mean-threshold)**2/(2*mean)).exp() for mean,threshold in budgets)
    risk = count_risk+300000000*((spread-64).exp()+(-q(190)).exp()+(spread-12000).exp())
    assert risk < q("6e-11")
    print("热端点自由能用于局部提议后的失败概率上界",risk.str(10),flush=True)

    phase = run_path(str(root/"ch39-first-arithmetic-plaquette-winding.py"))
    p = 23226
    phi = p*phase["angle"](11431)
    alpha = 2*q.pi()/3
    e = 2*(1-(phi-alpha).cos())/p
    reverse = 2*(1-(phi+alpha).cos())/p
    split = reverse-e
    assert 0 < e < reverse < q(".00013")
    assert 6-e > split/2
    # chi <= 1-eta 强制多付一个路径二次型下界。
    delta = q("1e-9")*split/2
    margin = 3*delta/((1+e/8)*(1+(e+delta)/8))
    orientation = margin/(2*heat_max)-spread-q(16).log()
    assert orientation > 4800
    print("目标扇区 chi <= 1-10^-9 的概率负对数下界",orientation.str(11))

    # 切面 |chi| <= 1/2 与两种等边取向分开，使用原始非负势给绝对下界。
    signed_area = 3*q(".7")*q(".98")**2-q(".1")*2*q("1.02")
    assert 2*signed_area/(q(3).sqrt()*3*q("1.02")**2) > q(".5")
    barrier = min(4*q(".02")**2,2*q(".1")**2)
    trial_scale = q("26.81676714125598")
    half_width, eps = q("1e-8"), q("1.005e-7")
    joint = balance["joint"]
    model = joint["flow"]["gauge"]["shortest"]["model"]
    center = model["target"](11431,trial_scale)
    allowance = eps+(q(".09")+2*center*(2*half_width).exp())*half_width
    contrast = max(q(joint["leads"]["terminal_modes"](k)[1]) for k in (1,2,3))
    witness = 3*reverse+3*contrast*allowance**2+p*q("2e-16")
    assert 1/(1-(-q(".5")).exp()) < 3
    bottleneck = (barrier/2-witness)/heat_max-spread-16-(16*3*5/half_width).log()
    assert bottleneck > q("2e13")
    print("负取向条件平衡的单步切面权比负对数下界",bottleneck.str(11))


if __name__ == "__main__":
    verify()
