"""在两档矢量状规范交点的整段验证轨迹上推进 Higgs 四次耦合。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify():
    running = run_path(str(Path(__file__).with_name("ch39-two-tier-dirac-running.py")))
    inputs = [(arb("125.25")/arb("246.22"))**2/2, arb(".165"), arb(".170")]
    values = inputs.copy()
    log_mass_factor = arb(0)
    first_possible = None
    first_certain = None
    def step(time,dt,gauge,y):
        nonlocal first_possible, first_certain, log_mass_factor
        pi = arb.pi()
        gy,g2 = 12*pi*gauge[0]/5, 4*pi*gauge[1]
        trace = 3*y[0]+3*y[1]+y[2]
        scale = 1/(16*pi*pi)
        linear = scale*(4*trace-9*g2-3*gy)
        const = scale*(arb(3)/8*(2*g2**2+(g2+gy)**2)-2*(3*y[0]**2+3*y[1]**2+y[2]**2))
        quadratic = 24*scale
        pad = arb(".04")*dt.upper()
        log_mass_factor += dt*scale*(12*(values[0]+arb(0,pad))+2*trace-arb(9)/2*g2-arb(3)/2*gy)
        def advance(endpoint):
            reach = endpoint + arb(0,pad)
            speed = quadratic*reach**2+linear*reach+const
            assert abs(dt*speed) < pad
            return endpoint+dt*speed
        if first_possible is None and not values[0]-pad > 0:
            first_possible = time.lower()
        for i,value in enumerate(values):
            # 对同一外部轨迹，标量 ODE 的解按初值严格保序。
            lower = advance(value.lower()).lower()
            upper = advance(value.upper()).upper()
            values[i] = lower.union(upper)
        if first_certain is None and values[0] < 0:
            first_certain = (time+dt).upper()
    _,mu,_ = running["verify"](step)
    assert first_possible is not None and first_certain is not None
    assert values[0]<0 and values[1]<0<values[2]
    mz = arb("91.1876")
    mass_factor = log_mass_factor.exp()
    mass_ratio = arb("125.25")**2*mass_factor/(2*mu**2)
    assert arb("-.061") < values[0] < arb("-.056")
    assert mz*first_possible.exp() > arb("8.0e5") and mz*first_certain.exp() < arb("1.14e6")
    print("二次质量参数乘数及其与交点平方之比",
          *[(v.lower().str(10),v.upper().str(10)) for v in (mass_factor,mass_ratio)])
    print("树级低能初值及交点四次耦合区间",inputs[0].str(11),
          values[0].lower().str(11),values[0].upper().str(11))
    print("首次零点能量界 / GeV",(mz*first_possible.exp()).str(12),(mz*first_certain.exp()).str(12))
    print("低端临界初值两侧的终点符号",values[1].str(8),values[2].str(8))
    print("零边界所需树级质量区间 / GeV",
          *[(arb("246.22")*(2*x).sqrt()).str(10) for x in inputs[1:]])


if __name__=="__main__":
    verify()
