"""十九点原生 Weyl 高模的磁性阈值差与必须另加的逆耦合分裂。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpz_mat

crossing = run_path(str(Path(__file__).with_name("ch39-unit-seventy-five-crossing.py")))
soft = crossing["soft"]


def verify():
    ctx.prec = 224
    h = acb.zeta_zero(11431).imag
    base, neutral, _ = soft["vector_pencil"](h)
    unweighted = base+neutral
    keep = [i for i in range(19) if i != soft["pole"]]
    pinned = fmpz_mat([[int(unweighted[i,j].real.unique_fmpz()) for j in keep] for i in keep])
    tau = int(pinned.det())
    assert tau == 41501049051635045807
    trace = sum(unweighted[i,i].real for i in range(19))
    assert trace == 270
    mean = trace/18
    decrement = arb(23)/25*arb(2)**-20/(h*h+soft["point"]["study"]["offset"])
    w = 1-decrement
    all_scales = (mean**18/(w**18*tau)).log()/2
    assert all_scales < arb("1.78633")

    # The pinned actual kernel is bounded by the pinned unit graph at
    # stiffness w. Compute a sharper bound for the prescribed shell scale.
    poles = [w*z.real for z,k in pinned.charpoly().complex_roots() for _ in range(k)]
    assert len(poles)==18 and all(v>0 for v in poles)
    radius = soft["geo"]["theta"](804)/soft["geo"]["theta"](12)
    fermion_spread, scalar_total = arb(0), arb(0)
    def fermion_cdf(y):
        return (3*y-y**3)/2
    for j in range(256):
        a,b = arb(j)/256,arb(j+1)/256
        z = radius**2*(1-a*a)/4
        pinned_sum = sum((1+z/v).log()/2 for v in poles)
        spread = pinned_sum-9*(1+z/mean).log()
        fermion_spread += (fermion_cdf(b)-fermion_cdf(a))*spread
        scalar_total += (b**3-a**3)*pinned_sum
    assert fermion_spread < arb("1.307") and scalar_total < arb("5.5")
    em,weak = arb("127.930"),arb(".23129")
    initial = 3*(1-weak)*em/5 - 3*weak*em + 2/arb(".1180")
    Lmax = (arb("1e17")/arb("91.1876")).log()
    baseline = initial+Lmax/(5*arb.pi())
    necessary = -(baseline+scalar_total/(5*arb.pi()))*arb.pi()*5/18
    residual = baseline+scalar_total/(5*arb.pi()) + arb(18)/(5*arb.pi())*fermion_spread
    assert residual < -arb("8.7")
    # For species-dependent overall mass-squared factors, the log derivative
    # of each of the 18 finite-subtraction thresholds is bounded by 1/2.
    stiffness_ratio = ((necessary-all_scales)/9).exp()
    assert stiffness_ratio > arb("2.20")
    # In C+J/105, the gauge and Higgs heavy terms have opposite, equal
    # coefficients. The scalar angular weight is dominated by the vector
    # weight since y^3 <= (8*y-y^3)/7 on the unit interval.
    j0 = 5*3*(1-weak)*em/5-3*weak*em-2/arb(".1180")
    Lmin = (arb("1e13")/arb("91.1876")).log()
    projected = initial+j0/105-Lmin/(105*arb.pi())
    cost = arb(27)/(7*arb.pi())
    assert (cost-(arb(128)/35+arb(1)/5)/arb.pi()).contains(0)
    worst = projected+cost*all_scales
    separation = -worst*35/106
    required_tree_defect = -projected/cost
    assert worst < -arb("8.8939") and separation > arb("2.936")
    assert required_tree_defect > arb("9.0303")
    print("Grounded tree cofactor, all-scale/Fermi-shape spread upper",tau,
          all_scales.upper().str(11),fermion_spread.upper().str(11))
    print("Scalar exposure upper, necessary Weyl spread, residual C upper",
          scalar_total.upper().str(11),necessary.lower().str(11),residual.upper().str(11))
    print("Necessary fermion stiffness ratio with fixed Higgs calibration",stiffness_ratio.lower().str(11))
    print("All common scales: projected residual upper, inverse spread lower, necessary tree defect",
          worst.upper().str(11),separation.lower().str(11),required_tree_defect.lower().str(11))
    return residual.upper()


if __name__ == "__main__":
    verify()
