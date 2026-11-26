"""三态环的全零点规范筛选和自由计次的锐代价区间。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

absolute = run_path(str(Path(__file__).with_name("ch39-fermion-absolute-matching.py")))



def vacuum(angle):
    hop = arb(1)/2000
    squared = [1+hop**2+2*hop*((angle+2*k*arb.pi())/3).cos() for k in range(3)]
    return -4*sum(x*x*x.log() for x in squared)/(16*arb.pi()**2)


def gauge_window():
    hop = arb(1)/2000
    gauge = absolute["gauge"]
    absolute["ground"]["flow_bounds"]()
    time = arb("26.81676714125598")
    tol = arb(2)**-21
    model = gauge["shortest"]["model"]
    ratios = gauge["shortest"]["ratios"]
    base = gauge["gauge_state"](gauge["word"],time,(18,6,3))
    mixing = [a+b for a,b in gauge["events"]["mixing"]["displacement"](27)]
    target = model["target"](11431,time)
    errors = [1/(x*v)-target+arb(0,m.upper()) for x,v,m in zip(base,ratios,mixing)]
    rates = [arb(lo).union(arb(hi)) for lo,hi in
             (("-.055","-.054"),(".036",".038"),(".089",".09"))]
    projection = [z*v for z,v in zip([5,3,-8],ratios)]
    assert sum(projection).contains(0)
    derivative = sum(z*y for z,y in zip(projection,rates))
    assert derivative < -arb(".84")
    shift_bound = 8*((1+hop**3)/(1-hop**3)).log()/(8*arb.pi()**2)
    width = (sum(abs(z)*(abs(e).upper()+tol) for z,e in zip(projection,errors))+shift_bound)/abs(derivative).lower()
    assert width < arb(".000012")
    domain = arb(0, width.upper())
    angles = arb(0).union(arb.pi())
    def shift(offset):
        return absolute["shifts"](arb(42),offset,time,angles,size=3,hop=hop)
    shifts = shift(domain)
    color = target+errors[2]+rates[2]*domain-shifts[2]
    assert color+tol < model["target"](11430,time+width)
    assert color-tol > model["target"](11432,time-width)
    offset = arb("4.3e-7")
    witness = [e+rate*offset+target*(1-(-2*offset).exp())-s/v for e,rate,s,v in
               zip(errors,rates,shift(offset),ratios)]
    worst = max(abs(z).upper() for z in witness)
    assert worst < tol
    print("全部编号的尺度半宽",width.str(11),"任意相位下 11431 三支残差上界",worst.str(11))


def cost_window():
    hop = arb(1)/2000
    theta = absolute["angle"](11431)
    alpha = 2*arb.pi()/3
    count = 23226
    delta = alpha-count*theta
    assert 0 < delta < theta/2
    best = vacuum(delta)
    upper = vacuum(delta+theta)-best
    rival = 92905
    lower = (best-vacuum(rival*theta-4*alpha))/(rival-count)
    cutoff = count+(best-vacuum(arb(0)))/lower
    assert 122935 < cutoff < 122936 and 122936*theta+delta < 6*alpha
    # 只有误差小于 delta 的后续整数才可能使代价下界提高。
    contenders = []
    for band in (1,2,4,5):
        lo, hi = (band*alpha-delta)/theta, (band*alpha+delta)/theta
        # 在浮点区间端点上向外加一格，包含恰相等的候选。
        for p in range(max(count+1,int(lo.lower().floor().unique_fmpz())-1),min(122935,int(hi.upper().ceil().unique_fmpz())+1)+1):
            distance = abs(p*theta-band*alpha)
            if p == rival:
                assert distance < delta
                contenders.append(p)
            else:
                assert distance > delta
    assert contenders == [rival]
    # 前半段每一步的弦差，由倒数第二步与起始端的正弦下界控制。
    slope_lower = 4*hop**3/(8*arb.pi()**2*(1+3*hop**2/4+3*hop**4+hop**6+2*hop**3))
    early_sine = min((delta+3*theta/2).sin().lower(),alpha.sin().lower())
    assert 2*slope_lower*(theta/2).sin()*early_sine > upper
    assert arb("2.16785e-26") < lower < arb("2.16787e-26")
    assert arb("4.06658e-20") < upper < arb("4.06659e-20")
    print("锐下端",lower.str(16),"锐上端",upper.str(16),"长计数截点",cutoff.str(14))
    return lower,upper


if __name__ == "__main__":
    ctx.prec = 180
    gauge_window()
    cost_window()
