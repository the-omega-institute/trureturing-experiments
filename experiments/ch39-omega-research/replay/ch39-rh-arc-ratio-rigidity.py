"""Enclose the mass ratio for two shell records, uniformly in positive realizations."""
from fractions import Fraction as Q
from math import prod
from pathlib import Path
from runpy import run_path
import json

directory = Path(__file__).parent
integer = run_path(str(directory/"ch39-rh-clock-integer-replay.py"))
Bounds, cell = integer["Bounds"], integer["read_cell"]


def envelope(a, count, ref, weights):
    t = a*a
    indices = tuple(range(1,count+1))
    b = prod((a-j)/(a+j) for j in indices)
    eta = 2*b/(1+b*b)
    coefficients = [2*j*j*eta/(t-j*j)*prod(Q(j+k,j-k) for k in indices if k!=j) for j in indices]
    # Independently check partial fractions at enough rational test points for
    # the cleared polynomial identity (degree at most 2*count+2 in u).
    for u in range(count+2,3*count+5):
        z = -Q(1,u*u)
        f = 1/(1+t*z)
        d = prod(Q(u-j,u+j) for j in indices)
        assert f-sum(c/(1+j*j*z) for j,c in zip(indices,coefficients))==f*eta*(d+1/d)/2
    radius = 2*abs(eta)+sum(map(abs,coefficients))/2**96
    base = sum((v*w/(v+Bounds(t)) for v,w in zip(ref,weights)),Bounds(0))
    return Bounds(base.lo-radius,base.hi+radius)


def main():
    data = json.loads((directory/"ch39-rh-clock-integer-data.json").read_text())["inputs"]
    ref = list(map(cell,data["reference"]))
    pivot = cell(data["pivot"][0])
    weights = [cell(x)/(Bounds(1)-pivot) for x in data["rho"]]
    spans = []
    for target, count, lo, hi in [(Q(33,34),38,"4.4335047516","4.4335047528"),
                                 (Q(14,17),41,"13.36660636879","13.36660636882")]:
        l, r = Q(lo), Q(hi)
        assert envelope(l,count,ref,weights).lo > target
        assert envelope(r,count,ref,weights).hi < target
        spans.append(Bounds(l,r))
    ratio = (spans[1]/spans[0])**2
    assert Q("9.089666741")<ratio.lo<=ratio.hi<Q("9.089666747")
    assert 6/(3-Q(2237,1000)) < ratio.lo  # sqrt(5) < 2.237
    assert Q(2237,1000)**2>5
    print("Every scalar or unit-vector matrix response requires squared ratio in (9.089666741,9.089666747).")


if __name__=="__main__":
    main()
