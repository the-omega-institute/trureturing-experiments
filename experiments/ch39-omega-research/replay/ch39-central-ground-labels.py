"""Optimize the central-input gauge shell and both clock directions at all labels."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb

central = run_path(str(Path(__file__).with_name("ch39-central-source-shell-test.py")))
floating = central["floating"]
leads = central["leads"]


def verify():
    ceiling = floating["ceiling"]
    cap = ceiling["verify"]()
    minimum, xs, t0, _ = central["central_landscape"]()
    pen = floating["penalty"]
    score = pen["competition"]["score"]
    free = pen["doublet"]["joint"]["contact"]["clock"]["free_ground"]
    ref = free(score(11431),2)
    advantage = ref-free(cap,2)
    modes = [tuple(map(arb,leads["terminal_modes"](count))) for count in (1,2,3)]
    allowance = arb(".000035")
    radius = arb(".15")
    spread = arb(".1435")*radius-abs(xs[0]-xs[2]).upper()
    target = floating["gauge"]["shortest"]["model"]["target"]
    mean = sum(xs)/3
    boundary = min((target(9199,t0+radius)-mean).lower(),(mean-target(14001,t0-radius)).lower())
    boundary -= arb(".09")*radius
    assert spread > 0 and boundary > 0
    assert min(d for m,d in modes)*spread**2/2 > allowance+advantage
    assert 3*min(m for m,d in modes)*boundary**2 > allowance+advantage
    screen = arb(".00740")
    assert arb(".007300684") < screen and free(screen,2)-ref > allowance
    pulse = ceiling["pulse"]
    candidates = {}
    # All prospective winners outside the two tested labels can be discarded
    # by transmission alone unless they clear this one-sided threshold.
    for start in range(9200,14001,32):
        for n,z in enumerate(acb.zeta_zeros(start,min(32,14001-start)),start):
            a,b = pulse["port_outputs"](pulse["generator"](z.imag),2,impulse=True)
            forward = abs((-7*a[0]+24*b[0])/25)**2
            if not forward < screen:
                assert forward > screen
                candidates[n] = free(forward,2)-ref
    print("central clock candidates",sorted(candidates),flush=True)
    for count,(m,d) in enumerate(modes,1):
        levels = {n:e+minimum(n,m,d) for n,e in candidates.items()}
        winner = min(levels,key=lambda n:float(levels[n].mid()))
        runner_floor = min([allowance]+[e.lower() for n,e in levels.items() if n!=winner])
        assert winner == (11431,10762,10762)[count-1]
        assert runner_floor-levels[winner] > arb(("2.52e-6","1.50e-6","2.07e-6")[count-1])
        # Every mixed-sign score is <= the forward score. A pure forward
        # orientation in the original order has zero singlet contact.
        print("central global winner",count,winner,"label margin >",
              (runner_floor-levels[winner]).lower().str(9),flush=True)


if __name__ == "__main__":
    verify()
