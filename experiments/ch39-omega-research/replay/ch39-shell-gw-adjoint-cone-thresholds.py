"""Adjoint band slopes, a variational density, and an eight-threshold cone on the two-clock torus."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb

first = run_path(str(Path(__file__).with_name("ch39-shell-gw-first-adjoint-ground.py")))
full, shell, q = first["full"], first["shell"], arb


def velocity():
    single,_,_,_ = shell["transported_spectra"]()
    eigenvalues = []
    for i,x in enumerate(np.linalg.eigvalsh(full["model"]["array"](single))):
        center = q(float(x))
        lo,hi = center-q("1e-13"),center+q("1e-13")
        assert full["count"](full["rows"](single),lo) == i
        assert full["count"](full["rows"](single),hi) == i+1
        eigenvalues.append(lo.union(hi))
    assert eigenvalues[1]-eigenvalues[0] > q("3.69")
    assert eigenvalues[-1]-eigenvalues[0] < 6
    eye = shell["yuk"]["identity"](6)
    projector = eye
    for other in eigenvalues[1:]:
        projector = projector*(single-other*eye)/(eigenvalues[0]-other)
    pts = shell["pts"]
    currents = []
    for source in range(6):
        a,b = pts[source].index(3),pts[source].index(4)
        if b != (a+1)%3:
            continue
        dest = next(i for i in range(6) if sum(v*w for v,w in zip(pts[i],pts[source])) == 25)
        currents.append((single[dest,source]*projector[source,dest]).imag)
    assert len(currents) == 3
    for current in currents:
        assert q("-.012977790") < current < q("-.012977788")
    # Equality of these oriented currents is exact by the cyclic spatial rotation.
    slope = 4*q(3).sqrt()*abs(currents[0])
    assert q(".08991275") < slope < q(".08991278")
    return slope


def eight_thresholds(slope):
    Q,E = first["gauge_target"]()
    tmax = q("1e-18")
    # On |p| <= t = |r|, the nine low levels move by <= 4t. The
    # orthogonal complement was >= lambda_s+3.69 and stays >3 above this band.
    assert q("3.69")-8*tmax > 3
    # Compression error <= 8t², Schur correction <= (4t)²/3.
    per_level = 8+q(16)/3
    delta = q(0,(2*per_level*tmax).upper())
    slopes = [slope*i/2+delta for i in (1,1,2,2,2,3,3,4)]
    def masses(mu):
        return [(125**2+mu**2*d).sqrt() for d in slopes]
    lo,hi = q("60410.25"),q("60410.28")
    assert sum((Q/m).log() for m in masses(lo)) > E
    assert sum((Q/m).log() for m in masses(hi)) < E
    group_ranges = [(masses(lo)[i].lower(),masses(hi)[i].upper()) for i in (0,2,5,7)]
    for (low,high),guess in zip(group_ranges,(12809,18114,22185,25617)):
        assert guess < low and high < guess+1
    assert lo*(3/tmax).sqrt() > Q
    # In the annular part t >= 1e-19 all 45 upper levels are sub-Planckian.
    assert 6+8*tmax < 7
    assert (125**2+hi**2*7/q("1e-19")).sqrt() < q("5.06e14")
    first["chiral_window"](q(0),q(0),tmax)
    # For the full cone (including angles arbitrarily close to zero) the second
    # copy is unnecessary for this lower bound and no UV cutoff is imposed.
    assert tmax**2/(2*q.pi()**2) > q("5.066e-38")
    assert (tmax**2-q("1e-19")**2)/(2*q.pi()**2) > q("5.015e-38")
    print("mu=M sqrt(|r|) in",(float(lo),float(hi)),"GeV; 8 active masses by group",
          [(float(l),float(u)) for l,u in group_ranges])


def main():
    slope = velocity()
    # One lobe pair at each of (0,0) and (pi,pi), related by (p+pi,pi-r).
    density = (slope/2)**2/(3*q.pi()**2)
    assert density > q("6.8259e-5")
    assert density/4 > q("1.7064e-5")
    print("Lowest cusp slope",slope.str(12),"adjoint inversion density above",density.lower().str(10))
    eight_thresholds(slope)


if __name__ == "__main__":
    main()
