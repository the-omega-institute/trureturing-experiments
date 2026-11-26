"""Find the fewest records when each species and its reference use one orbit."""
from bisect import bisect_left
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb

source = run_path(str(Path(__file__).with_name("ch39-relative-chord-word-source.py")))
events, short, q = source["events"], source["short"], source["q"]
counts = (169, 312, 138, 6)
reference = (18, 6, 3)
orbits = ((36,1,1), (31,2,0), (13,5,2), (14,12,12), (29,7,2),
          (9,2,0), (37,1,0), (12,12,11), (15,12,7))


def reference_candidates(weights, low_bounds=events["low"], high_bounds=events["high"]):
    ordered = sorted(weights.values())
    viable = list(weights)
    for name in ("e", "u", "d", "s", "c", "mu", "b", "tau", "t"):
        low, high = low_bounds[name], high_bounds[name]
        rest = []
        for key in viable:
            mass = weights[key]
            minimum = (mass*low.numerator+low.denominator-1)//low.denominator
            idx = bisect_left(ordered, minimum)
            if idx < len(ordered) and ordered[idx]*high.denominator <= mass*high.numerator:
                rest.append(key)
        viable = rest
    return set(viable)


def verify():
    rows = [[v*v for v in word] for word in source["chord_words"]]
    mass = source["weighted_ports"](rows, counts)
    assert all(sum(v) == 5**17 for v in rows)
    assert len(set((*orbits, reference))) == 10
    for copies in range(1, 39):
        redraw = events["orbit_weights_from_ports"](copies, mass)
        fixed = [events["orbit_weights_from_ports"](copies, v) for v in rows]
        batch = {k:sum(n*data[k] for n,data in zip(counts, fixed, strict=True)) for k in redraw}
        assert sum(batch.values()) == 625*5**(17*copies)
        for name, weights in (("redraw", redraw), ("batch", batch)):
            feasible = reference_candidates(weights)
            if copies < 38:
                assert not feasible, (name, copies, feasible)
                continue
            assert reference in feasible
            ys = [Fraction(weights[k], weights[reference]) for k in orbits]
            assert all(events["low"][n] < v < events["high"][n]
                       for n,v in zip(events["names"], ys, strict=True))
            state = [arb(v)/mass[0] for v in mass[1:]]+[q(v) for v in ys]
            errors = short["residuals"](state, q("26.8166617100"))
            assert max(errors) < q("3.18e-7")
            ratios = [(q(v).sqrt()/q(events["centers"][k]))
                      for k,v in zip(events["names"], ys, strict=True)]
            print(name, "38-copy reference choices", len(feasible), "Yukawa ratios",
                  min(ratios).lower().str(10), max(ratios).upper().str(10),
                  "all-mixing gauge residual bounds", *(v.upper().str(10) for v in errors))


if __name__ == "__main__":
    verify()
