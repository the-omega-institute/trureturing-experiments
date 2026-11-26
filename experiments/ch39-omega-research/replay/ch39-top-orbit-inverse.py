"""Determine the top orbit and the two possible references in the book input box."""
from bisect import bisect_left, bisect_right
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
single = run_path(str(root/"ch39-single-orbit-flavor-source.py"))
source = single["source"]
shell, events, q = source["shell"], single["events"], single["q"]
solver = shell["shortest"]
lead = run_path(str(root/"ch39-unit-lead-gauge-shell.py"))
top = (13, 5, 2)
time = q("26.8166617100")


def residuals_at(gauge, logs, scale):
    state = [*gauge, *(x.exp() for x in logs)]
    for _ in range(128):
        state = solver["reference"]["advance"](state, scale/128, 14)
    target = shell["model"]["target"](11431, scale)
    return [1/(r*v)-target for r,v in zip(solver["ratios"], state)]


def orbit_choices(weights, reference):
    return {name:sorted(k for k,v in weights.items()
                        if events["low"][name]*weights[reference] <= v <= events["high"][name]*weights[reference])
            for name in events["names"]}


def second_ground(gauge, ys, speeds, displacement, scale):
    target = shell["model"]["target"]
    arithmetic = arb(0, q("1.3e-16"))
    uncertainty = [arb(0,(a+q(".003")*b).upper()) for a,b in displacement]
    mean_speed = sum(speeds)/3
    for side in (-1, 1):
        t = scale+side*q("4e-9")
        raw = residuals_at(gauge, ys, t)
        errors = [x+e+arithmetic for x,e in zip(raw, uncertainty, strict=True)]
        # Positive residual slopes hold throughout this bracket; it suffices
        # to check both endpoints for a uniform tolerance bound.
        target_at_t = target(11431,t)+arithmetic
        assert all(v+2*target_at_t > 6 for v in speeds)
        assert all(abs(e) < arb(2)**-21 for e in errors)
        mean = sum(errors)/3
        contrast = sum((r-mean)*(v-mean_speed) for r,v in zip(errors,speeds,strict=True))
        for i in (1,2,3):
            m,d = lead["terminal_modes"](i)
            gradient = 3*m*mean*(mean_speed+2*target_at_t)+d*contrast
            assert side*gradient > 0


def verify():
    counts = single["counts"]
    ports = [[v*v for v in w] for w in source["chord_words"]]
    mass = source["weighted_ports"](ports,counts)
    gauge = [arb(v)/mass[0] for v in mass[1:]]
    redraw = events["orbit_weights_from_ports"](38,mass)
    tables = [events["orbit_weights_from_ports"](38,row) for row in ports]
    batch = {k:sum(n*t[k] for n,t in zip(counts,tables,strict=True)) for k in redraw}
    _,_,logs,cons,beta,_ = shell["scale_constraints"](time)
    displacement = events["mixing"]["displacement"](27)
    speeds = shell["short"]["floating"]["flow_bounds"]()
    for name,weights in (("redraw",redraw),("batch",batch)):
        ordered = sorted(weights.values())
        references = {r for r,w in weights.items() if all(
            bisect_left(ordered, w*events["low"][k]) < bisect_right(ordered, w*events["high"][k])
            for k in events["names"] if k != "t")}
        assert references == {(17,3,2),(18,3,2),(19,5,2),(18,6,3)}
        ref = single["reference"]
        reference_logs = [q(Fraction(weights[k],weights[ref])).log() for k in single["orbits"]]
        errors = residuals_at(gauge, reference_logs, time)
        # Use the weak and colour channels. Gauge inputs are identical, so
        # only derivatives with respect to the nine log-squared Yukawas enter.
        row,_,margin = cons[1]
        weak = errors[1]-beta[1]*errors[2]
        offset = weak+sum(row[3+j]*(logs[j]-reference_logs[j]) for j in range(9) if j != 2)
        assert not row[5].contains(0)
        band = ((reference_logs[2]-(offset+arb(0,margin.upper()))/row[5])/2).exp()
        assert q(".85594") < band < q(".85811")
        viable = []
        for r in sorted(references):
            choices = orbit_choices(weights,r)
            for k,w in weights.items():
                value = q(Fraction(w,weights[r])).sqrt()
                if not value.overlaps(band):
                    continue
                if (r,k) == ((19,5,2),(16,8,2)):
                    # Retain all feasible choices for the other eight inputs.
                    selected = [l[0] for l in choices.values()]
                    selected[2] = k
                    anchor = [q(Fraction(weights[z],weights[r])).log() for z in selected]
                    bound = errors_at_fixed_top(gauge, weights, r, choices, anchor, row, beta[1])
                    assert abs(bound) > margin
                else:
                    viable.append((r,k))
        assert viable == [((17,3,2),top),((18,6,3),top)]
        r = (17,3,2)
        choices = orbit_choices(weights,r)
        selected = [l[0] for l in choices.values()]
        selected[2] = top
        assert len(set((*selected,r))) == 10
        ys = [q(Fraction(weights[z],weights[r])).log() for z in selected]
        t = q("26.8166614942" if name == "redraw" else "26.8166614927")
        second_ground(gauge,ys,speeds,displacement,t)
        print(name, "top orbit", top, "couplings at the two references",
              *(q(Fraction(weights[top],weights[r])).sqrt().str(12) for r,_ in viable))


def errors_at_fixed_top(gauge, weights, ref, choices, anchor, row, beta):
    errors = residuals_at(gauge, anchor, time)
    bound = errors[1]-beta*errors[2]
    for j,(name,opts) in enumerate(choices.items()):
        if name != "t":
            possible = [q(Fraction(weights[k],weights[ref])).log() for k in opts]
            bound += row[j+3]*(min(possible).union(max(possible))-anchor[j])
    return bound


if __name__ == "__main__":
    verify()
