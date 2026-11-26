"""Bound the second singular value throughout the first two scalar modes."""
from itertools import combinations
from pathlib import Path
from runpy import run_path
from flint import arb, acb_mat

m = run_path(str(Path(__file__).with_name("ch39-shell-gw-asymmetric-bands.py")))
m["cuts"][2] = ".015882"
b,c = m["orb"]["orthonormalize"](m["images"](2))
size = list(combinations(range(3),2))


def exterior(x,y):
    return acb_mat([[x[i,k]*y[j,l]-x[i,l]*y[j,k] for k,l in size] for i,j in size])


if __name__ == "__main__":
    cols = [exterior(b,b),(exterior(b,c)+exterior(c,b))/arb(2).sqrt(),exterior(c,c)]
    gram = acb_mat([[m["hs"](x,y) for y in cols] for x in cols])
    threshold = arb(".28")
    assert m["m"]["count"](m["m"]["rows"](gram),threshold**2) == 0
    # HS-normalized matrices satisfy sum of squared products of singular values
    # <= (2r²+r⁴)/(1+r²)² when r = s₂/s₁, by increasing s₃ up to s₂
    # in the numerator and dropping it from the denominator.
    r = arb(".2")
    assert (2*r*r+r**4)/(1+r*r)**2 < threshold**2
    print("two-mode second/largest singular-value ratio exceeds",r.str(4))
