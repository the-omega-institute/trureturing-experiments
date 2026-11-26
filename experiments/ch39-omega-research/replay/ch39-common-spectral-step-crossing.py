"""第 11431 枚的共同单调两档质量函数与规范耦合有限减法一圈交点。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def crossing(fermion_copies=1):
    model = run_path(str(Path(__file__).with_name("ch39-native-center-parity-crossing.py")))
    ctx.prec = 256
    spectra, _ = model["setup"]()
    counts = {0:0, 1:2, 4:1, 2:0, 3:0, 6:1}
    for q,row in spectra.items():
        assert len(row) == 18
        assert sum(int(v < arb("7.66")) for v in row) == counts[q]
        assert sum(int(v > arb("7.72")) for v in row) == 18-counts[q]

    pi = arb.pi()
    fermion = lambda m, t: model["sheet"](m,t,False)
    vector = lambda m, t: model["sheet"](m,t,True)
    scalar = model["scalar"]
    em, weak = arb("127.930"), arb(".23129")
    initial = [(1-weak)*em*3/5, weak*em, 1/arb(".1180")]
    low = [3*fermion_copies*sum(counts[q]*row[i] for q,row in zip(model["charges"],model["weyl"]))
           for i in range(3)]
    assert all((v-w).contains(0) for v,w in zip(low,[fermion_copies*arb(16)/5,6*fermion_copies,5*fermion_copies]))

    def exposure(kernel, mass_squared, log_scale):
        # 三类角积分具有单位角权。log(1+t)<=t 对整个 M_Z 减项给出包围。
        bound = (-2*log_scale).exp()/(8*mass_squared)
        return kernel(mass_squared,arb(1)) - arb(0).union(bound.upper())

    def high_baseline(log_scale):
        h,v,f = [18*exposure(kernel,arb(1),log_scale) for kernel in (scalar,vector,fermion)]
        return [initial[i]+(-model["beta"][i]*log_scale + arb(7)*n*v/2 -
                 model["higgs"][i]*h - 3*fermion_copies*sum(row[i] for row in model["weyl"])*f)/(2*pi)
                for i,n in enumerate((0,2,3))]

    def crossing_data(log_scale):
        baseline = high_baseline(log_scale)
        delta = 2*pi*(baseline[1]-baseline[2])/(low[1]-low[2])
        adjusted = [a-v*delta/(2*pi) for a,v in zip(baseline,low)]
        return adjusted[0]-adjusted[2], delta, adjusted[2]

    def bracket(fn,left,right,steps):
        lo,hi = arb(left),arb(right)
        sign = 1 if fn(hi)>0 else -1
        assert sign*fn(lo)<0 and sign*fn(hi)>0
        for _ in range(steps):
            mid=(lo+hi)/2
            value=sign*fn(mid)
            if value>0:
                hi=mid
            else:
                assert value<0
                lo=mid
        return lo.union(hi)

    L = bracket(lambda t:crossing_data(t)[0],29,30,42)
    _, delta, common = crossing_data(L)
    def fermion_residual(x):
        return exposure(fermion,(-2*x).exp(),L)-exposure(fermion,arb(1),L)-delta
    x = bracket(fermion_residual,7,16,35)
    mu = model["mz"]*L.exp()
    light_mass = mu*(-x).exp()
    ratio = 18/(16+2*(-2*x).exp())
    assert arb("5e14")<mu<arb("7e14") and arb("29")<common<arb("32")
    assert arb("1.124999") < ratio < arb("1.125")
    print("每个重费米表示的 Weyl 计数", fermion_copies)
    print("交点 GeV、共同逆耦合",mu.str(15),common.str(14))
    print("较低档质量 GeV、对数质量比、十八层平方质量和比",
          light_mass.str(14),x.str(14),ratio.str(17))


if __name__ == "__main__":
    for copies in (1,2):
        crossing(copies)
