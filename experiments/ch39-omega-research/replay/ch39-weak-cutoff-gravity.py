"""A finite proper-time cutoff bounds induced gravity for every nonnegative scalar curvature coupling."""
from pathlib import Path
import runpy
from flint import arb

m = runpy.run_path(str(Path(__file__).with_name("ch39-curvature-spectrum-rigidity.py")))
q,A = m["q"],m["weyl"]
pi = arb.pi()
scale = 1000*q("27.0455").exp()
mp = q("1.220890e19")
needed = (mp/scale)**2
# Count each tower by rank on the graph. Vectors plus the one zero-mode
# compensator have a nonpositive Einstein heat coefficient in this matter model.
ns = 2*10 + 20
nf = 2*8 + 16
positive = q(ns+nf)/6
def budget(x):
    return ((2*A/3)*x**3+19*x**2)/(2*pi)

def check():
    x = q("5.19350")
    assert ns == 40 and nf == 32
    fraction = positive*budget(x)/needed
    # xi >= -b gives scalar coefficient at most 1/6+b per real direction.
    b = (needed/budget(x)-positive)/ns
    copies = 6*needed/budget(x)
    lo,hi = q(1),q(1000)
    for _ in range(50):
        z=(lo+hi)/2
        if positive*budget(z) < needed:
            lo=z
        else:
            hi=z
    assert positive*budget(lo)<needed
    print("proper-time fraction at weak cutoff",fraction.upper().str(14))
    print("necessary maximum negative scalar xi magnitude",b.lower().str(14))
    print("required total scalar or Weyl rank if xi>=0",copies.lower().str(14))
    print("minimum cutoff / Q at original ranks",lo.lower().str(14))

if __name__ == "__main__":
    check()
