"""One-loop SU(2) continuation above the fitted crossing with the stated weak representations."""
from pathlib import Path
import runpy
from flint import arb

short = runpy.run_path(str(Path(__file__).with_name("ch39-short-arc-vacuum.py")))
q, r, sh = (short[k] for k in ("q","r","shell"))
pi = arb.pi()
levels = []
for _, _, d, _ in sh["edges"]:
    length = r * (sh["rho"] * sh["theta"](d)).sqrt() * short["fsbox"][short["band"](d)]
    fundamental = pi / length
    assert fundamental > 1
    n = 1
    while (n * fundamental).lower() < 10:
        levels.append(n * fundamental)
        n += 1
# Each individual scalar edge level has multiplicity two in S_R.
lows = sorted(m.lower() for m in levels)
highs = sorted(m.upper() for m in levels)
assert len(levels) > 9

def accumulated(x):
    # Integrate b2=-1/2+N_R/2. Rank-38 interlacing gives the upper envelope.
    # N_R >= max(18,N_D) gives the lower envelope: replace the first nine
    # doubled Dirichlet levels by the eighteen levels already active at Q.
    lower = q("8.5")*x.log() + sum((x/m).log() for m in highs[9:] if m < x)
    upper = q("18.5")*x.log() + sum((x/m).log() for m in lows if m < x)
    return lower.lower().union(upper.upper())

def crossing(alpha_inverse):
    targets = 2*pi*(q("40.6142").union(q("40.9906"))-alpha_inverse)
    ranges = []
    for side in ("early","late"):
        lo, hi = q(1), q(10)
        assert accumulated(hi) > targets
        for _ in range(38):
            mid = (hi+lo)/2
            bounds = accumulated(mid)
            if ((bounds.upper() < targets.lower()) if side=="early"
                else not (bounds.lower() > targets.upper())):
                lo=mid
            else:
                hi=mid
        ranges.append(lo if side=="early" else hi)
    return ranges[0].lower().union(ranges[1].upper())

if __name__ == "__main__":
    for label, value in (("alpha_2=1",1),("inverse alpha_2=0",0)):
        x = crossing(value)
        print(label,"energy / Q",x.lower().str(13),x.upper().str(13))
