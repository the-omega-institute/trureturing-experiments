"""Realize the top Yukawa by orbit-invariant pooling of repeated quaternary port records."""
from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent
s = run_path(str(root/"ch39-shortest-matching-rotor.py"))
mixing = run_path(str(root/"ch39-ckm-gauge-displacement.py"))
q = s["q"]


def orbit_weights(copies):
    weights = [v*v for v in s["match_coord"]]
    orbits = defaultdict(int)
    for outcome in product(range(4),repeat=copies):
        key = tuple(sorted([outcome.count(j) for j in (1,2,3)],reverse=True))
        orbits[key] += prod(weights[j] for j in outcome)
    assert sum(orbits.values()) == (5**21)**copies
    return [orbits[key] for key in sorted(orbits)]


def event_ratio(weights,assignment):
    numerator = sum(w for w,a in zip(weights,assignment) if a == 1)
    denominator = sum(w for w,a in zip(weights,assignment) if a == 2)
    if numerator and denominator:
        return Fraction(numerator,denominator)


def verify():
    ctx.prec = 128
    center = Fraction(".8616")
    enclosure = s["Enclosure"]()
    # For a rule to work uniformly in mixing it must work in the diagonal case.
    # The preceding shortest-word certificate bounds the only possible label
    # and the low-energy top Yukawa at this depth.
    low,high = [(center*Fraction(f))**2 for f in (".9921",".9947")]
    for copies in (1,2,3):
        weights = orbit_weights(copies)
        selected = set()
        for assignment in product((0,1,2),repeat=len(weights)):
            value = event_ratio(weights,assignment)
            if value is not None and low <= value <= high:
                selected.add(value)
        for value in selected:
            top = q(value).log()
            jet = enclosure.central_jet(s["match_coord"],top)
            assert not enclosure.possible_jet(s["match_label"],jet,arb(0))
        print(copies,"copies, candidate ratios in the necessary top-Yukawa window",len(selected))
    assignment = (0,2,1,1,2,2,0,1,1,0,1)
    squared = event_ratio(orbit_weights(4),assignment)
    factor = q(squared).sqrt()/q(center)
    assert q(".99364757") < factor < q(".99364758")
    a,*imag = s["match_coord"]
    ys = [z*z for z in s["model"]["yukawa_centers"]]
    ys[2] = q(squared)
    time = q("26.81621480455")
    state = [(arb(i)/a)**2 for i in imag]+ys
    for _ in range(128):
        state = s["reference"]["advance"](state,time/128,14)
    change = [a+b for a,b in mixing["displacement"](27)]
    target = s["model"]["target"](s["match_label"],time)
    errors = [abs(1/(u*r)-target)+d for u,r,d in zip(state,s["ratios"],change)]
    assert all(err<s["tolerance"] for err in errors)
    print("Four-copy rate ratio",squared,"top Yukawa",q(squared).sqrt().str(13),
          "all-mixing residual <",max(err.upper() for err in errors).str(9))


if __name__ == "__main__":
    verify()
