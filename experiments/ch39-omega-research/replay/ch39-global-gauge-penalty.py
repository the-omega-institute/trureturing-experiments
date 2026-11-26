"""Compare a mean gauge constraint with a positive three-channel energy on every zero."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb

root = Path(__file__).parent
doublet = run_path(str(root/"ch39-unordered-port-doublet.py"))
competition = run_path(str(root/"ch39-two-tick-global-competition.py"))
gauge = doublet["joint"]["contact"]["gauge"]


def verify():
    maximum, endpoint_mass, acceptance = doublet["verify"]()
    free = doublet["joint"]["contact"]["clock"]["free_ground"]
    bare = free(maximum,2)
    time = arb("26.81676714125598")
    target = gauge["shortest"]["model"]["target"]
    errors = gauge["gauge_match"](gauge["word"],11431,time,(18,6,3))
    allowance = arb("1.005e-7")
    assert all(e < allowance for e in errors)
    center = target(11431,time)
    lower, upper = target(2994,time),target(14840,time)
    low_probability = (center-upper)/(lower-upper)
    assert lower > center > upper and arb(".0394094") < low_probability < arb(".0394095")
    scores = [competition["score"](n) for n in (2994,14840)]
    energy_advantage = bare-low_probability*free(scores[0],2)-(1-low_probability)*free(scores[1],2)
    assert energy_advantage > arb(".0000087")
    variance = (lower-center)*(center-upper)
    assert variance>35

    neighbor = min((target(11430,time)-center).lower(),(center-target(11432,time)).lower())
    assert neighbor>2*allowance
    outside = min((target(10814,time)-center).lower(),(center-target(12093,time)).lower())
    tail_cost = 3*(outside-allowance)**2
    ground_ceiling = bare+3*allowance**2
    assert tail_cost-ground_ceiling > arb(".0417")
    window_gap = free(maximum-arb(".00002295"),2)-bare
    separation = min(window_gap.lower(),(tail_cost-ground_ceiling).lower())
    assert separation > arb("4e-5")
    # The RH reciprocal-zero sum bounds every finite reference normalization.
    normalization = 1+arb.const_euler()/2-(4*arb.pi()).log()/2
    height = acb.zeta_zero(11431).imag
    weight_lower = 1/((arb(1)/4+height**2)*normalization)
    bright = endpoint_mass*weight_lower*acceptance/(8*625)
    gap = bright*separation/(separation+acceptance/2)
    assert gap > arb("1.86e-15")
    print("exterior mixture",low_probability.str(10),"energy gain",energy_advantage.str(9),
          "gauge variance",variance.str(9))
    print("outside energy margin", (tail_cost-ground_ceiling).lower().str(9),
          "gap for every finite cutoff containing 11431 >",gap.lower().str(9))


if __name__ == "__main__":
    verify()
