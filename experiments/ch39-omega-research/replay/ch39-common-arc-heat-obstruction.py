"""Uniform two-moment obstruction for a right tower on the same dynamic arcs and endpoint space."""
from pathlib import Path
import runpy

import numpy as np
from flint import arb

root = Path(__file__).resolve().parent
overlap = runpy.run_path(str(root/"ch39-scalar-overlap-bound.py"))
source = overlap["source"]
q = overlap["q"]
r = q(400)

def heat_constants():
    # Bounds for osc(g), sup|x*g'(x)| and sup|g|, with g=x*f'(x).
    # For x>=24, x^3 exp(-x) is decreasing.
    tail = 4*q(24)**3*(-q(24)).exp()
    assert tail<q(".00001")
    granges = [arb(0,tail),arb(0,tail)]
    slopes = [tail,tail]
    for i in range(24*128):
        x = q(i).union(q(i+1))/128
        z = (-x).exp()
        for j,(g,s) in enumerate(((-x*z,x*(x-1)*z),(x*(1-x)*z,x*(1-3*x+x*x)*z))):
            granges[j] = granges[j].union(g)
            slopes[j] = max(slopes[j],abs(s).upper())
    shell = source["shell"]
    mass = sum(4*shell["rho"]*shell["theta"](d) for _,_,d,_ in shell["edges"])
    mesh = q(".01")
    hs2 = mass*(q("1/3")+mesh**2/6)
    endpoint_trace = mass*mesh/2
    errors = []
    for g,s in zip(granges,slopes):
        c = 2*max((g.upper()-g.lower()).upper(),s)
        errors.append(2*c*hs2+abs(g).upper()*endpoint_trace)
    return errors

def all_bounds(H):
    values = np.linalg.eigvalsh(np.array([[complex(H[i,j].mid()) for j in range(38)] for i in range(38)]))
    result = []
    for j,v in enumerate(values):
        a,b = arb(str(v))-q(".002"),arb(str(v))+q(".002")
        assert a>0 and source["inertia"](H,a)<=j<source["inertia"](H,b)
        result.append(a.union(b))
    return result

def exclude(low,weights,errors):
    T,W = overlap["ranges"](low)
    lower = all_bounds(source["static"]([w.lower() for w in weights]))
    upper = all_bounds(source["static"]([w.upper() for w in weights]))
    ys = [lo.lower().union(hi.upper())/r**2 for lo,hi in zip(lower,upper)]
    # The common Dirichlet sums cancel. Both vertex heat readings pay the error twice.
    S = sum((-v).exp() for v in ys)+arb(0,2*errors[0].upper())
    M = sum(v*(-v).exp() for v in ys)+arb(0,2*errors[1].upper())
    e = q(1).exp()
    assert S-M>17 and S+(e-1)*M<19 and 0<M<1/e
    n = 18
    loss = n*(1-(-(-2*W.upper()/n).exp()).exp()-(-2*T.lower()).exp())
    high = S.lower()-n+loss
    assert 0<high<1/e
    high_m = -high*high.log()
    assert M.upper()-high_m < q(".3")*(-q(".3")).exp()
    # Each low y is now below .3, so y exp(-y) >= .3/(exp(.3)-1)*(1-exp(-y)).
    minimum = high_m + q(".3")/(q(".3").exp()-1)*loss
    assert minimum > M.upper()
    print("right static S,M",S.lower().str(12),S.upper().str(12),M.lower().str(12),M.upper().str(12))
    print("necessary high heat and total moment, contradiction margin",
          high.lower().str(12),minimum.lower().str(12),(minimum-M.upper()).lower().str(12))

if __name__=="__main__":
    errors = heat_constants()
    print("one graph: errors in heat and first moment",*(v.upper().str(12) for v in errors))
    exclude(overlap["profile"]["low"],[arb(0)]*13,errors)
    exclude(source["mass_bounds"](),source["source_bounds"](),errors)
