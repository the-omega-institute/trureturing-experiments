"""Test resolved gauge-matched ground labels over the entire original input box."""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb

root = Path(__file__).parent
leads = run_path(str(root/"ch39-unit-lead-gauge-shell.py"))
floating = leads["floating"]
window = run_path(str(root/"ch39-floating-window-refinement.py"))


def verify():
    tolerance = arb(2)**-21
    terminal_error = arb("1.3e-16")
    widened = tolerance+terminal_error
    outer, _ = window["bands"](window["scale_bounds"](widened),widened)
    labels = window["base"]["certified_labels"](*outer)
    assert labels == (10815,12092)
    rates = floating["flow_bounds"]()
    target = floating["gauge"]["shortest"]["model"]["target"]
    test = arb(27)
    center = target(11431,test)
    # The ratio of two height targets is independent of the matching scale.
    shift = max(abs((center/target(n,test)).log()/2).upper() for n in labels)
    assert shift < arb(".04871")
    anchor = window["q"](window["fixed"]["L"])
    assert 23 < anchor-arb(1)/8-shift and anchor+arb(1)/8+shift < 33
    modes = [tuple(map(arb,leads["terminal_modes"](k))) for k in (1,2,3)]
    mean_term = 3*max(m for m,d in modes)*(widened+terminal_error+shift*abs(sum(rates)/3))**2
    contrast_term = max(d for m,d in modes)*sum(
        (2*widened+shift*abs(ri-rj))**2 for ri,rj in combinations(rates,2))/3
    cost = mean_term+contrast_term
    assert cost < arb(".000027")
    doublet = floating["penalty"]["doublet"]
    maximum, _, _ = doublet["verify"]()
    free = doublet["joint"]["contact"]["clock"]["free_ground"]
    gap = free(maximum-arb(".00002295"),2)-free(maximum,2)
    assert gap-cost > arb(".000013")
    print("all sources: matched competing label is beaten by 11431 by at least",
          (gap-cost).lower().str(9),"trial cost below",cost.upper().str(9))


def flavor():
    fit = run_path(str(root/"ch39-radial-gain-two-angle-fit.py"))
    mixture = run_path(str(root/"ch39-clock-flavor-mixture.py"))
    model = fit["model"]
    gain = arb(18)/25
    centers = model["shell"]["centers"]
    widths = [model["shell"]["sine_radius"]]*3+[model["shell"]["phase_radius"]]
    for n in (11431,10762):
        response = model["fixed"]["native"]["plane_response"](acb.zeta_zero(n).imag)
        def values(params):
            matrix = model["observables"]((gain,*params),response=response)
            return mixture["readouts"](mixture["direct_observables"](matrix))
        common = values([arb("2.561167"),arb("3.08245")])
        errors = [abs(v-c)/w for v,c,w in zip(common,centers,widths,strict=True)]
        assert max(e.upper() for e in errors) < arb(".84")
        print(n,"same two flavor parameters",*(v.str(9) for v in common))
        if n==10762:
            original = values([arb("2.5605152"),arb("3.0832980")])
            assert original[0] < centers[0]-widths[0]
            outer,solution = fit["isolate"](("2.56","3.08"),response,gain)
            fit["prove_least"](response,gain,outer)
            held = values(solution)
            assert all(abs(v-c)<w for v,c,w in zip(held,centers,widths,strict=True))
            print("10762 least-contact calibration: held-out sine and phase",held[2].str(9),held[3].str(9))


if __name__ == "__main__":
    flavor()
    verify()
