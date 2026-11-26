"""Test all shortest common-rotor gauge preparations across the full energy window."""
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

root = Path(__file__).parent
words = run_path(str(root/"ch39-hurwitz-gauge-words.py"))
volume = run_path(str(root/"ch39-volume-contenders.py"))
model, q = volume["gauge"],volume["q"]


def verify():
    ctx.prec = 128
    _,_,logs = volume["input_bounds"]()
    ratios = model["ratios"]
    margins = []
    for quat, letters in words["hits"](16):
        a,*imag = quat
        inv = [(arb(a)/v)**2 for v in imag]
        gaps = []
        def observe(step,x,rate):
            if step >= 23*64:
                speeds = [v/r for v,r in zip(rate,ratios)]
                assert speeds[0] < speeds[1] < speeds[2]
        model["running"]["integrate"]([abs(arb(v)/a) for v in imag],33*64,observe)
        # For each fixed trajectory the first-minus-other channel gaps decrease.
        # A positive gap to its right and a negative gap to its left exclude every L.
        for offset in ("-.021","+.023","+.003","-.008"):
            time = volume["anchor"] + q(offset)
            end,_,_ = model["evolve"](inv,logs,time,8192,diagonal_quarks=True)
            g1,g2,g3 = [v/r for v,r in zip(end,ratios)]
            diffs = [g1-g2,g1-g3]
            positive = max(v.lower() for v in diffs)
            negative = max((-v).lower() for v in diffs)
            gaps.append(min(positive,negative))
        best = max(gaps)
        assert best/2 > q(".00017224")
        margins.append(best)
    assert len(margins) == 4
    print("Shortest words, global lower bound for half-spread of normalized couplings",
          (min(margins)/2).str(10))


if __name__ == "__main__":
    verify()
