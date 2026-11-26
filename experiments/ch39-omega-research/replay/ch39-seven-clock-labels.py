"""在两个弱极点上达到禁止直接四次读数的最小味标签数。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import ctx

source = run_path(str(Path(__file__).with_name("ch39-discrete-clock-charge.py")))


def verify():
    ctx.prec = 180
    power = 23226
    step, extra = divmod(power,5)
    assert (step,extra) == (4645,1)
    powers, weights = source["chain"](step)
    powers.append(step+1)
    assert 4*max(powers) < power
    # 四份末级乘积与一份再乘边界值的乘积，最后一次有一个独立误差。
    total_error_weight = 5*sum(weights[1:])+1
    assert total_error_weight == power-5
    states = [("A",0),("A",step),("A",2*step),
              ("B",2*step),("B",3*step),("B",4*step),("B",0)]
    source["check_ring"](states, [step]*4+[step+1], len(powers)-1, total_error_weight)
    multiplicities = [sum(site==k for site,charge in states) for k in ("A","B")]
    beta = [Fraction(-16,3)+n*Fraction(4,3) for n in multiplicities]
    assert beta == [Fraction(-4,3),0]
    # 临界站点含六个复双基本标量和八个 Dirac 基本双态。
    adjoint, fundamental = 2, Fraction(3,4)
    pure = -Fraction(34,3)*adjoint**2
    fermions = (Fraction(20,3)*adjoint+4*fundamental)*4
    scalars = (Fraction(2,3)*adjoint+4*fundamental)*6
    assert pure+fermions+scalars == 46
    print("极点一圈系数",beta,"临界站点两圈纯自身规范系数",pure+fermions+scalars)


if __name__ == "__main__":
    verify()
