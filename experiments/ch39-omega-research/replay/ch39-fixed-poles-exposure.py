"""Dependence of the scalar threshold on a positive weight interpolation at fixed poles."""
from fractions import Fraction as Q
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).resolve().parent
a = run_path(str(root/"ch39-cutoff-thresholds.py"))
B = a["Bounds"]


def ball(x):
    def convert(r):
        return arb(r.numerator)/r.denominator
    return convert(x.lo).union(convert(x.hi))


def compute():
    ctx.prec = 480
    records = a["data"]("ch39-cutoff-threshold-inverses.json")
    base, alt = [a["exact_fit"](records[i]) for i in (2,1)]
    ids = sorted(set(base)|set(alt))
    assert len(ids)==49 and set(base)^set(alt)=={29,455}
    poles = list(map(ball,a["kinetic"](ids)))
    h = a["data"]("ch39-tick-separators.json")["prefix"][0][0]
    mean = poles[0]*ball(B(Q(2**800,(h+1)**2),Q(2**800,h*h)))
    knots = [-mean*j*j for j in a["OFFSETS"]]
    # Common denominator of the difference of the two Stieltjes functions.
    # The numerator vanishes at 44 known arguments, leaving degree at most 4.
    centers = [45,81,127,210,473]
    nodes = [poles[ids.index(i)] for i in centers]
    vals = []
    for i,c in zip(centers,nodes):
        dw = ball(base.get(i,B(0))-alt.get(i,B(0)))
        numerator = dw*prod(1-c/d for k,d in zip(ids,poles) if k!=i)
        vals.append(numerator/prod(c-t for t in knots))
    coefficients = [v/prod(c-d for j,d in enumerate(nodes) if j!=i)
                    for i,(v,c) in enumerate(zip(vals,nodes))]
    def difference(x):
        poly = sum(v*prod(x-d for j,d in enumerate(nodes) if j!=i)
                   for i,v in enumerate(coefficients))
        return poly*prod(x-t for t in knots)/prod(1-x/d for d in poles)
    p,modes,q02,_ = a["physical_inputs"]()
    p,q02 = ball(p),ball(q02)
    modes = list(map(ball,modes))
    ends = [[ball(w.get(i,B(0))) for i in ids] for w in (alt,base)]
    mixed = [[(1-t)*l+t*r for l,r in zip(*ends)] for t in (arb(0),arb(3)/4)]
    rhoends = [[(1-p)*u for u in w] for w in mixed]
    rho = [x.union(y) for x,y in zip(*rhoends)]
    bracket = (104358280169550,104358280169551)
    for t in (Q(0),Q(3,4)):
        weights = {i: B(t)*base.get(i,B(0))+B(1-t)*alt.get(i,B(0)) for i in ids}
        for m,sign in zip(bracket,(1,-1)):
            assert a["endpoint"](m,ids,weights,passive_above=True)[:2] == (sign,144)
    # The 48-site endpoint B omits coordinate 455. In the 49-site family
    # it is counted in the kinetic mean for every positive interpolation.
    assert set(ids)-set(alt)=={455}
    assert bracket[0]-a["MASS_BRACKETS"][1][1] > 1046813876750
    mr = arb(bracket[0]).union(arb(bracket[1]))
    higgs = (125/mr)**2
    q = q02/mr**2
    count = sum(c<q for c in poles)
    assert count==3 and poles[count]>q
    slopes = []
    mass_counts = arb(25*count)
    for lam,n in zip(modes,a["MODES"]):
        def f(x):
            values = [x-lam*p-higgs-lam*x*sum(w/(x-c) for w,c in zip(row,poles)) for row in rhoends]
            return values[0].union(values[1])
        for j in range(count+1):
            # Locate uniform left and right bounds separately, allowing any t.
            start = poles[j-1].upper() if j else arb(0)
            end = poles[j].lower()
            l,r=start,end
            for _ in range(82):
                m=(l+r)/2
                if f(m)<0:
                    l=m
                else:
                    r=m
            lo=l
            l,r=start,end
            for _ in range(82):
                m=(l+r)/2
                if f(m)>0:
                    r=m
                else:
                    l=m
            hi=r
            assert f(lo)<0 and f(hi)>0
            x=lo.union(hi)
            assert x>start and x<end
            if j==count:
                assert x>q
                continue
            assert x<q and x*mr**2>1000**2
            dx = 1+lam*sum(w*c/(x-c)**2 for w,c in zip(rho,poles))
            measure = 11*n*lam*(1-p)*difference(x)/(x*dx)
            slopes.append(measure)
            mass_counts += n*(1-higgs/(x*dx))
    tangent = mr*sum(slopes)/(22*mass_counts)
    assert arb("3.44e-29")<tangent and tangent<arb("3.45e-29")
    print("Parent mass slope in GeV per unit mixing",tangent.str(16))
    print("Exposure derivative",sum(slopes).str(18))
    print("by gap", [sum(slopes[j::count]).str(10) for j in range(count)])
    return sum(slopes)


if __name__ == "__main__":
    compute()
