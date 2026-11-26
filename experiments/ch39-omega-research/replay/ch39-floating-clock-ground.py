"""Jointly minimize over the shell label, the state and the classical RG scale."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb

root = Path(__file__).parent
penalty = run_path(str(root/"ch39-global-gauge-penalty.py"))
ceiling = run_path(str(root/"ch39-two-phase-pulse-ceiling.py"))
gauge = penalty["gauge"]


def flow_bounds():
    mixing = gauge["events"]["mixing"]
    model = mixing["model"]
    running, ratios, loop = model["running"], model["ratios"], model["loop"]
    low, high, logs = mixing["volume"]["input_bounds"]()
    box = [a.union(b) for a,b in zip(low,high)]
    initial = gauge["rotor_inputs"](gauge["word"], (18,6,3))
    assert all(x.contains(1/u) for x,u in zip(box, initial[:3], strict=True))
    assert all(y.contains(u.log()) for y,u in zip(logs, initial[3:], strict=True))
    speed_bands = [arb(lo).union(arb(hi)) for lo,hi in
                   (("-.055","-.054"),(".036",".038"),(".089",".09"))]
    count = 0
    def observe(step, g, y, traces):
        nonlocal count
        count += 1
        if count <= 23*128:
            return
        rates = [-2*(beta+sum(b*u-d*t for b,d,u,t in zip(br,dr,g,traces))/loop)/(loop*r)
                 for beta,br,dr,r in zip(running["BETA"],running["B"],running["D"],ratios)]
        assert rates[2]-rates[0] > arb(".1435") and all(abs(v) < arb(".09") for v in rates)
        assert all(box.contains(v) for box,v in zip(speed_bands,rates,strict=True))
        # The true squared gauge and Yukawa velocities enter one more derivative
        # of the inverse gauge flow. A convex quark mixing average is <= its trace.
        dg = [abs(v*r*u*u) for v,r,u in zip(rates,ratios,g)]
        total = 3*(traces[0]+traces[1])+traces[2]
        trace_speeds = []
        for family,row in enumerate(running["C"]):
            other = traces[1-family] if family<2 else 0
            gauge_sum = sum(c*u for c,u in zip(row,g))
            trace_speeds.append(sum(u*(3*(u+other)+2*(total+gauge_sum))/loop
                                    for u in y[3*family:3*family+3]))
        for br,dr,r in zip(running["B"], running["D"],ratios):
            acc = 2*sum(abs(b)*u+abs(d)*v for b,d,u,v in zip(br,dr,dg,trace_speeds))/(loop**2*r)
            assert acc < 1
    model["evolve"](box, logs, arb(33), 33*128, observe_flavors=observe)
    assert count == 33*128
    return speed_bands


def verify():
    score_limit = ceiling["verify"]()
    maximum, endpoint_mass, acceptance = penalty["doublet"]["verify"]()
    free = penalty["doublet"]["joint"]["contact"]["clock"]["free_ground"]
    advantage = free(maximum,2)-free(score_limit,2)
    assert 0 < advantage < arb(".000011193")
    speed_bands = flow_bounds()
    time = arb("26.81676714125598")
    allowance = arb("1.005e-7")
    errors = gauge["gauge_match"](gauge["word"],11431,time,(18,6,3))
    assert all(e < allowance for e in errors)
    target = gauge["shortest"]["model"]["target"]
    center = target(11431,time)
    incumbent = 3*allowance**2

    def outside_distance(radius, below, above):
        drift = allowance + arb(".09")*radius
        return min((target(below,time+radius)-center-drift).lower(),
                   (center-drift-target(above,time-radius)).lower())

    cutoff = arb(".04695")
    exterior = outside_distance(cutoff,10814,12093)
    far_spread = arb(".1435")*cutoff-2*allowance
    assert exterior > 0 and far_spread > 0
    outside_cost = min((3*exterior**2).lower(),(far_spread**2/2).lower())
    margin = (outside_cost-incumbent)/2-advantage
    assert margin > arb("1e-7")

    # If any label attains a gauge cost <= incumbent, the gauge spread forces
    # its scale to lie in this small neighborhood of the trial scale.
    narrow = (2*incumbent.sqrt()+2*allowance)/arb(".1435")
    neighbor = outside_distance(narrow,11430,11432)
    assert narrow < arb(".000004") and neighbor > incumbent.sqrt()

    # The incumbent label has a unique minimizer. In its localized interval
    # all three residuals are increasing, and the sum of squares is strictly convex.
    low, high = center*(-2*narrow).exp(), center*(2*narrow).exp()
    slope = 2*low-arb(".09")
    residual = allowance+arb(".09")*narrow+high-low
    curvature = 6*(slope**2-residual*(1+4*high))
    scale_radius = allowance/slope
    assert curvature > 250 and scale_radius < arb("1.55e-8")
    # At the optimal classical scale the gauge cost of any other label is
    # greater. The outside-window labels even have a macroscopic cost here.
    optimum_outside = outside_distance(scale_radius,10814,12093)
    assert optimum_outside > 0
    window_gap = free(maximum-arb(".00002295"),2)-free(maximum,2)
    exterior_gap = (3*optimum_outside**2-incumbent)/2-advantage
    assert exterior_gap > window_gap > arb(".00004")
    normalization = 1+arb.const_euler()/2-(4*arb.pi()).log()/2
    height = acb.zeta_zero(11431).imag
    bright = endpoint_mass*acceptance/(8*625*(arb(1)/4+height**2)*normalization)
    gap = bright*window_gap/(window_gap+acceptance/2)
    assert gap > arb("1.86e-15")
    # A pure forward state is dark for the singlet at every label.
    witness_shift = (target(14840,time)/center).log()/2+arb(".00082")
    assert 23 < time+witness_shift < 33
    mismatch = center-target(14840,time+witness_shift)
    witness_cost = sum((allowance+abs(mismatch+v*witness_shift))**2 for v in speed_bands)
    witness_gain = free(maximum,2)-free(penalty["competition"]["score"](14840),2)
    assert witness_gain > arb(".0158")*witness_cost
    print("14840 defeats all 11431 scales for strengths through 0.0158")
    print("bare energy advantage below",advantage.upper().str(9),
          "exterior margin for strength >= 1/2 above",margin.lower().str(9),flush=True)
    print("unique classical scale within",scale_radius.upper().str(9),
          "curvature above",curvature.lower().str(9),"quantum gap above",gap.lower().str(9))


if __name__ == "__main__":
    verify()
