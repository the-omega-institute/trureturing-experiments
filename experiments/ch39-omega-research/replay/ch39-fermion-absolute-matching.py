"""给原 128 字母输入加入完整费米子多重态的一圈有限阈值。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent
ground = run_path(str(root/"ch39-floating-clock-ground.py"))
gauge = ground["gauge"]
angle = run_path(str(root/"ch39-first-arithmetic-plaquette-winding.py"))["angle"]


def shifts(ratio, time_offset, center_time, delta, *, size=5, hop=None):
    # 质量以 M 为标尺，ratio=M/Q_*，低能减法点为 1 TeV。
    if hop is None:
        hop = arb(1)/2
    def log_sum(phase, momentum_ratio):
        total = arb(0)
        for j in range(size):
            mass_squared = 1+hop**2+2*hop*((phase+2*j*arb.pi())/size).cos()
            z = momentum_ratio**2/mass_squared
            linear = z/5
            total += (linear-3*z**2/140).union(linear)
        return total
    sums = [(log_sum(phase,time_offset.exp()/ratio)
             -log_sum(phase,(-center_time).exp()/ratio))
            for phase in (delta,arb(0))]
    return [(b*sums[0]+c*sums[1])/(16*arb.pi()**2)
            for b,c in [(arb(4)/5,arb(8)/15),(arb(4)/3,0),(0,arb(4)/3)]]


def verify():
    ctx.prec = 160
    time = arb("26.81676714125598")
    tol = arb(2)**-21
    target = gauge["shortest"]["model"]["target"](11431,time)
    ratios = gauge["shortest"]["ratios"]
    state = gauge["gauge_state"](gauge["word"],time,(18,6,3))
    mixing = [a+b for a,b in gauge["events"]["mixing"]["displacement"](27)]
    residuals = [1/(v*x)-target+arb(0,m.upper()) for v,x,m in zip(ratios,state,mixing)]
    rates = [arb(lo).union(arb(hi)) for lo,hi in
             (("-.055","-.054"),(".036",".038"),(".089",".09"))]
    ground["flow_bounds"]()
    # 监测量在整个 23 <= L <= 33 时消去新增阈值。
    monitor = [5,-3,-2]
    weights = [m*v for m,v in zip(monitor,ratios)]
    raw_rate = sum(w*r for w,r in zip(weights,rates))
    assert raw_rate < -arb(".55") and sum(weights).overlaps(-arb(2)/5)
    budget = sum(abs(w)*(abs(x).upper()+tol) for w,x in zip(weights,residuals))
    wide_radius = (budget/arb(".55")).upper()
    fast_rate = -raw_rate.upper()+4*target*( -2*wide_radius).exp()/5
    radius = (budget/fast_rate).upper()
    assert radius < arb("1.8e-6")
    common_rate = 2*target*((-2*radius).exp().union((2*radius).exp()))
    speeds = [v+common_rate for v in rates]
    delta = abs(23226*angle(11431)-2*arb.pi()/3)
    # 原始逆耦合间的阈值差界，单位为 1/g^2。
    difference = -((delta/2).cos().log())/(6*arb.pi()**2)
    slope = speeds[0].upper()/speeds[1].lower()
    max_shift = (residuals[0].upper()+tol
                 -slope*(residuals[1].lower()-tol)+difference/ratios[0])/(1/ratios[0]-slope/ratios[1])
    weak_at_33 = shifts(arb(33),-radius,time,delta)[1]
    # 负时间增量单独用第一支即可排除。
    assert weak_at_33 > ratios[0]*(residuals[0].upper()+tol)+difference
    assert weak_at_33 > max_shift
    def residual_box(mass_ratio, offset):
        return [x+speed*offset-y/v for x,speed,y,v in
                zip(residuals,speeds,shifts(arb(mass_ratio),offset,time,delta),ratios)]
    # 选取一个固定时间移量，使全部原 CKM 误差窗均通过。
    raw = shifts(arb(42),arb(0),time,delta)[1].mid()
    offset = (raw*(1/ratios[0]+1/ratios[1])/(speeds[0].mid()+speeds[1].mid())).mid()
    assert abs(offset)<radius
    boxes = residual_box(42,offset)
    assert all(abs(x)<tol for x in boxes)
    print("匹配尺度半宽",radius.str(11),"第二支原始逆耦合移量上界",max_shift.str(11),
          "M/Q*=33 时下界",weak_at_33.lower().str(11),"42 时尺度移量",offset.str(11),
          "三支最大残差",max(abs(x).upper() for x in boxes).str(11))


if __name__ == "__main__":
    verify()
