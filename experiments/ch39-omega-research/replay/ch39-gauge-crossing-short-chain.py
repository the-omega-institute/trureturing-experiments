"""在两支规范量交叉处从未知零点编号准备短链，不用指定初始编号。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
dynamics = run_path(str(root/"ch39-local-winding-metastability.py"))
joint = dynamics["joint"]
q = arb


def verify():
    ctx.prec = 144
    nominal = q("26.81676714125598")
    eps = q("1.005e-7")
    gauge = joint["flow"]["gauge"]
    model = gauge["shortest"]["model"]
    target = model["target"]
    n = 11431
    height = acb.zeta_zero(n).imag
    radius = 2*eps/q(".1435")
    scale = nominal+q(0,radius.upper())
    center = target(n,scale)
    constant = model["fixed"]["OFFSET"]
    offset = q(constant.numerator)/constant.denominator
    assert 0 < offset < 31 and acb.zeta_zero(1).imag>14
    error = eps+radius*(q(".09")+2*center.upper())
    unit = center*(height**2+offset)
    bound = 3*height*(1+offset/14**2)**2/unit
    temperature = q("2e-17")
    prior_slope = temperature*(1+offset/14**2)/center
    low_gap = (target(n-1,nominal)-target(n,nominal))*(-2*(scale-nominal)).exp()
    high_gap = (target(n,nominal)-target(n+1,nominal))*(-2*(scale-nominal)).exp()
    assert high_gap>2*error
    for leads in (1,2,3):
        common,_ = [q(v) for v in joint["leads"]["terminal_modes"](leads)]
        left = 3*common*(low_gap-2*error)-bound
        right = 3*common*(high_gap-2*error)
        assert left>prior_slope
        gap = min((low_gap*left).lower(),(high_gap*right).lower())
        assert gap>q("4e-8")
        print("leads",leads,"root scale radius",radius.str(9),"short-chain label gap",
              gap.str(10),flush=True)
    total_weights=(1+q.const_euler()/2-(4*q.pi()).log()/2)*(q(1)/4+height**2)
    assert (q("4e-8")/temperature-total_weights.log())>q("1.9e9")
    assert q(100000)/2-23225>25000
    assert 2*q(25000)**2/100000>12000


if __name__=="__main__":
    verify()
