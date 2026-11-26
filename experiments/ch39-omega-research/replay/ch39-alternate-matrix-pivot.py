"""在相同标尺与矩阵动力学下比较两个已经登记的低能输入。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent
alternate = run_path(str(root/"ch39-cross-zero-cp-obstruction.py"))
currents = run_path(str(root/"ch39-two-conserved-current-pivot.py"))
scalar = run_path(str(root/"ch39-free-complex-messengers.py"))


def verify():
    ctx.prec = 160
    label, count = 11618,23555
    center = "26.802773689811964"
    word = alternate["word"]
    currents["verify"](label=label,center=center,word=word,discarded=None)
    gauge = alternate["gauge"]
    data = gauge["rotor_inputs"](word)
    mixing = gauge["events"]["mixing"]
    lo,hi,logs = mixing["volume"]["input_bounds"]()
    assert all(x.union(y).contains(1/v) for x,y,v in zip(lo,hi,data[:3],strict=True))
    assert all(b.contains(x.log()) for b,x in zip(logs,data[3:],strict=True))
    time = arb(center)
    mass_ratio = 42*(arb("26.81676714125598")-time).exp()
    target = gauge["shortest"]["model"]["target"](label,time)
    ratios = gauge["shortest"]["ratios"]
    state = gauge["gauge_state"](word,time)
    residual = [1/(v*x)-target+arb(0,(a+b).upper())
                for v,x,(a,b) in zip(ratios,state,mixing["displacement"](27))]
    def central_gauge(r, extra):
        offset = arb("4.2e-7")
        mass = (1-2*r).union(1+2*r)
        def kernel(momentum):
            z = (momentum/(mass_ratio*mass))**2
            return (z/5-3*z**2/140).lower().union((z/5).upper())
        common = (kernel(offset.exp())-kernel((-time).exp()))/(4*arb.pi()**2)
        speeds = [arb(s).union(arb(t)) for s,t in
                  (("-.055","-.054"),(".036",".038"),(".089",".09"))]
        errors = [abs(x+v*offset+target*(1-(-2*offset).exp())-common/ratio).upper()
                  for x,v,ratio in zip(residual,speeds,ratios)]
        assert max(errors)<arb(2)**-21
        print("Alternate central three-current residual",max(errors).str(10))
    scalar["verify"]("0",label=label,count=count,gauge_check=central_gauge)
    # 一般矩阵降到标量的粗预算只须重算新次数、零点对应的单位试探和振子上界。
    R = arb("3e-18")*(arb(25000)/42)**4
    r = arb("4e-5")
    tau = 6*r+9*r**2
    C = 4/(8*arb.pi()**2)
    rem = 3*C*tau**3/(6*(1-tau))
    delta = 2*arb.pi()/3 - count*scalar["free"]["absolute"]["angle"](label)
    budget = (2*rem + (C*r**2)**2/R + R*12*(delta/2).sin()**2/count
              + R*arb("5e-18")*3*arb.pi()*(count-1)/2)
    assert (budget/R).sqrt() < arb(".001381")
    print("Alternate full-matrix minimum",label,count,"positive orientation")


if __name__ == "__main__":
    verify()
