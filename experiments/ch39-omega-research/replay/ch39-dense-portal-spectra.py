"""Bound the individual active levels for 0 < d <= 1/2 in the zero-35 angle window."""
from pathlib import Path
from fractions import Fraction
import runpy
from typing import NamedTuple

from flint import arb

root = Path(__file__).resolve().parent
source = runpy.run_path(str(root / "ch39-heavy-inertia-quartic.py"))
count = source["eigen_count"]
q = source["q"]


class DenseSpectrum(NamedTuple):
    scale_lower: arb
    scale_upper: arb
    count_lower: int
    count_upper: int
    log_upper: list
    response_lower: arb
    mass_ratios: list


def spectra():
    h = q("1e-6")
    beta = arb.pi()**2*h*h/6
    floor = q(".05849")
    lower = [floor] * 308
    upper = [None] * 308
    for j in range(446):
        e = q(".05501") + q(Fraction(j, 200))
        bound_lo = e/(1+beta*(e+20))
        bound_hi = e/((1-h*h*e)**2-9*arb.pi()**2*h*(80+q(".1")*e))
        n_lo = count(e, q(0), 1)
        n_hi = count(e, q(".5"), -1)
        for k in range(n_hi, len(lower)):
            lower[k] = max(lower[k], bound_lo.lower())
        for k in range(n_lo):
            if upper[k] is None:
                upper[k] = bound_hi.upper()

    def exposure(r, xs):
        return sum(max(q(0), (r / x.sqrt()).log().upper()) for x in xs if x is not None)

    grid = [q("1.30")+q(Fraction(j,1000)) for j in range(151)]
    rlo = max(r for r in grid if exposure(r, lower) < q("48.23"))
    # A lower enclosure for exposure at given ratio.
    def lower_exposure(r):
        return sum(max(q(0), (r / x.sqrt()).log().lower()) for x in upper if x is not None)
    rhi = min(r for r in grid if lower_exposure(r) > q("49.02"))
    nlo = sum(x is not None and x<rlo**2 for x in upper)
    nhi = len(lower)-sum(x>rhi**2 for x in lower)
    ratios = [(lower[i]/rhi**2).lower() for i in range(nhi)]
    # At fixed i, a high bound on scalar log exposure and a low on y*(1-log y).
    logs = [-x.log()/2 for x in ratios]
    min_response = sum(x*(1-x.log()) for x in ratios[:nlo])
    boxes = [ratios[i].union(min(q(1), (upper[i]/rlo**2).upper())) for i in range(nlo)]
    return DenseSpectrum(rlo, rhi, nlo, nhi, logs, min_response, boxes)


if __name__ == "__main__":
    lo, hi, nlo, nhi, logs, response, _ = spectra()
    print("scale ratio",lo,hi,"counts",nlo,nhi,"common response",response.lower().str(12))
    print("six largest log bounds", [v.upper().str(7) for v in logs[:6]])
