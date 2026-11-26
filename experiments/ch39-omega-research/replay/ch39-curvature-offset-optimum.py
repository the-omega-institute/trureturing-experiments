"""Bound the first positive-curvature instability for every common mass offset y M^2."""
from pathlib import Path
import runpy
from flint import arb

pairs = runpy.run_path(str(Path(__file__).with_name("ch39-gapped-curvature-pairs.py")))
q, r, graph = (pairs[k] for k in ("q", "r", "graph"))
pi, rotations, degree = (pairs[k] for k in ("pi", "rotations", "degree"))
point = graph["point"]
j0 = q(".00030629")
Q = (1000*q("26.8220").exp()).union(1000*q("27.0455").exp())
mp = q("1.220890e19")
need = pi*(mp/Q)**2 - 2*q(".05").union(q("1.99779"))/3 - q(".000048797").union(j0)/4

def scan():
    start = q(1)/8
    grid = [start]
    step = q(2)**(q(1)/32)
    for _ in range(32*33):
        grid.append(grid[-1]*step)
    # For large y, h(y)=y J(y) <= 608 degree^3/(r^2 y).
    assert grid[-1] > 2*degree
    tail = (608*degree**3/(r**2*grid[-1]**2)).upper()
    upper = [tail]
    lower = [q(0)]
    for lo, hi in reversed(list(zip(grid,grid[1:]))):
        delta = (hi-lo)/r**2
        upper.append((upper[-1]+delta*point(lo,True,rotations).upper()).upper())
        lower.append((lower[-1]+delta*point(hi,False,rotations).lower()).lower())
    upper.reverse()
    lower.reverse()
    # Optimize h/(1-J/(4 need)) to keep the scalar's 1/4 exactly.
    corr = (1/(1-j0/(4*need.lower()))).upper()
    samples = [(lower[i]*y).lower() for i,y in enumerate(grid)]
    best = max(range(len(samples)),key=lambda i:samples[i])
    witness = samples[best]
    caps = [(hi*upper[i]*corr).upper() for i,hi in enumerate(grid[1:])]
    outside = max((start*j0*corr).upper(),(grid[-1]*tail*corr).upper())
    assert outside < witness
    active = [i for i,val in enumerate(caps) if not val < witness]
    maximum = witness.union(max(caps))
    arg = grid[active[0]].lower().union(grid[active[-1]+1].upper())
    base = 3*32   # (1/8)*2^3 = 1
    assert grid[base].contains(1)
    gain = (witness/(upper[base]*corr)).lower().union((max(caps)/lower[base]).upper())
    scale = ((Q.lower()/r)**2/need.upper()).lower().union(((Q.upper()/r)**2/need.lower()).upper())
    optimum = (maximum * scale).sqrt()
    energy_bound = (32*pi**2*Q.upper()**2*q("9.37987e-12")/(need.lower()-j0/4)).sqrt()
    assert optimum < energy_bound
    for label, val in (("test offset y",grid[best]),("maximizer y",arg),("optimized y J correction",maximum),
                       ("gain in critical R over y=1",gain),("optimal sqrt(R) GeV",optimum)):
        print(label,val.lower().str(14),val.upper().str(14))
    print("moment bound on sqrt(R) GeV",energy_bound.upper().str(14))

if __name__ == "__main__":
    scan()
