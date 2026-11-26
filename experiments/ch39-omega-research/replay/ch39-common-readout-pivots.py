"""两条不同的四元数词共用完全相同的 33 份端口分类和能量模型。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent
alternate = run_path(str(root/"ch39-alternate-matrix-pivot.py"))
gauge = alternate["alternate"]["gauge"]
events = gauge["events"]
removed, added = (18, 6, 3), (14, 10, 1)


def common_inputs(word):
    coords = (1, 0, 0, 0)
    letters = tuple(map(int, word))
    assert len(letters) == 128 and all((a ^ 1) != b for a,b in zip(letters, letters[1:]))
    for letter in letters:
        coords = gauge["words"]["multiply"](coords, gauge["words"]["GENERATORS"][letter])
    assert sum(x*x for x in coords) == 5**128
    weights = events["orbit_weights"](33,coords)
    masses,used = events["event_masses"](weights)
    table = (root/"ch39-common-yukawa-events.tsv").read_text().splitlines()
    assert any(line.split() == [*map(str,removed),"ref"] for line in table)
    assert added not in used
    ref = masses["ref"] - weights[removed] + weights[added]
    ys = {n:Fraction(masses[n],ref) for n in events["names"]}
    assert all(events["low"][n]<y<events["high"][n] for n,y in ys.items())
    gs = [abs(Fraction(x,coords[0])) for x in coords[1:]]
    assert all(lo<x<hi for x,(lo,hi) in zip(gs,gauge["words"]["BANDS"],strict=True))
    q = gauge["shortest"]["q"]
    return [q(x*x) for x in gs] + [q(y) for y in ys.values()]


def verify():
    ctx.prec = 160
    cases = [(11431, "26.81676714125598", "6.60e-7", gauge["word"]),
             (11618, "26.802773689811964", "1.68e-7", alternate["alternate"]["word"])]
    for label, log_scale, offset_text, word in cases:
        state = common_inputs(word)
        mixing = events["mixing"]
        lo,hi,logs = mixing["volume"]["input_bounds"]()
        assert all(x.union(y).contains(1/v) for x,y,v in zip(lo,hi,state[:3],strict=True))
        assert all(b.contains(x.log()) for b,x in zip(logs,state[3:],strict=True))
        center = arb(log_scale)
        for _ in range(128):
            state = gauge["shortest"]["reference"]["advance"](state,center/128,14)
        alternate["currents"]["verify"](label=label,center=log_scale,evolved_state=state)
        target = gauge["shortest"]["model"]["target"](label,center)
        ratios = gauge["shortest"]["ratios"]
        residual = [1/(v*x)-target+arb(0,(a+b).upper())
                    for v,x,(a,b) in zip(ratios,state,mixing["displacement"](27))]
        offset = arb(offset_text)
        mass_ratio = 42*(arb("26.81676714125598")-center).exp()
        r = arb("4e-5")
        mass = (1-2*r).union(1+2*r)
        def kernel(momentum):
            z = (momentum/(mass_ratio*mass))**2
            return (z/5-3*z**2/140).lower().union((z/5).upper())
        shift = (kernel(offset.exp())-kernel((-center).exp()))/(4*arb.pi()**2)
        rates = [arb(l).union(arb(h)) for l,h in (("-.055","-.054"),(".036",".038"),(".089",".09"))]
        errors = [abs(x+v*offset+target*(1-(-2*offset).exp())-shift/g).upper()
                  for x,v,g in zip(residual,rates,ratios)]
        assert max(errors) < arb(2)**-21
        print("Shared table, unique label",label,"full central residual below",max(errors).str(10))


if __name__ == "__main__":
    verify()
