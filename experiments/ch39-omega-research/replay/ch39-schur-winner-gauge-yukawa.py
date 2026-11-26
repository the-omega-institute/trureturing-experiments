"""Certify a common quaternion input at the minimum-gain shell zero."""
from fractions import Fraction
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent
events = run_path(str(root/"ch39-common-yukawa-port-copies.py"))
shortest = events["model"]
words = shortest["words"]
word = ("0404242552511431313024220521442030212131342034420352511443530433"
        "0522055220241342402535111240550351421314034155224000404253400213")


def rotor_inputs(sequence, discarded=None):
    ctx.prec = 144
    letters = list(map(int, sequence))
    assert len(letters) == 128 and all((a^1) != b for a,b in zip(letters,letters[1:]))
    coords = (1,0,0,0)
    for index in letters:
        coords = words["multiply"](coords,words["GENERATORS"][index])
    assert sum(x*x for x in coords) == 5**128
    for signs in product((-1,1), repeat=3):
        endpoint = (coords[0], *(v*s for v,s in zip(coords[1:],signs)))
        assert len(words["factor"](endpoint,128)) == 128
    for b,(lo,hi) in zip(coords[1:],words["BANDS"]):
        assert lo < abs(Fraction(b,coords[0])) < hi
    for copies in range(1,33):
        assert not events["possible_light_events"](events["orbit_weights"](copies,coords))
    weights = events["orbit_weights"](33,coords)
    masses,_ = events["event_masses"](weights)
    if discarded is not None:
        table = (root/"ch39-common-yukawa-events.tsv").read_text().splitlines()
        assert any(line.split() == [*map(str,discarded),"ref"] for line in table)
        masses["ref"] -= weights[discarded]
    ys = [Fraction(masses[n],masses["ref"]) for n in events["names"]]
    for n,z in zip(events["names"],ys):
        assert events["low"][n] < z < events["high"][n]
    return [(arb(v)/coords[0])**2 for v in coords[1:]] + [shortest["q"](v) for v in ys]


def gauge_state(sequence, time, discarded=None):
    state = rotor_inputs(sequence, discarded)
    for _ in range(128):
        state = shortest["reference"]["advance"](state,time/128,14)
    return state


def gauge_match(sequence, label, log_scale, discarded=None):
    time = arb(log_scale)
    state = gauge_state(sequence,time,discarded)
    target = shortest["model"]["target"]
    values = [1/(r*u) for r,u in zip(shortest["ratios"],state)]
    deltas = [a+b for a,b in events["mixing"]["displacement"](27)]
    errors = [abs(v-target(label,time))+d for v,d in zip(values,deltas)]
    for other in (label-1,label+1):
        assert all(abs(v-target(other,time))-d > arb(2)**-21 for v,d in zip(values,deltas))
    print("128-letter rotor, 33-port table,",label,"all-CKM residuals",
          *(e.upper().str(8) for e in errors))
    return errors


if __name__ == "__main__":
    errors = gauge_match(word,11431,"26.81676714125598",(18,6,3))
    assert all(e < arb("1.005e-7") for e in errors)
