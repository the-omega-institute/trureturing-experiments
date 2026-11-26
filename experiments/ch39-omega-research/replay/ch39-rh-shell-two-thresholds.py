"""Exact rational enclosures for reading two nine-point shell arcs at mass ratio sqrt(3)."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path
from runpy import run_path

directory = Path(__file__).parent
integer = run_path(str(directory/"ch39-rh-clock-integer-replay.py"))
Bounds, cell = integer["Bounds"], integer["read_cell"]


def main():
    points = [v for v in product(range(5), repeat=3) if sum(i*i for i in v)==17]
    arcs = sorted({Q(17+sum(i*j for i,j in zip(v,w)),34) for v,w in combinations(points,2)})
    assert len(points)==9 and len(arcs)==7 and max(arcs)==Q(33,34)
    a, b = Q(33,34), Q(14,17)
    assert a/(3-2*a)==Q(11,12) and a/(3-2*a)-b==Q(19,204)
    assert a*(1-b)/(b*(1-a))==Q(99,14)
    # Limiting uniform error across arbitrary positive probability mixtures.
    def polynomial(e):
        return 1156*e**2+2142*e-57
    assert polynomial(Q("0.02623907"))<0<polynomial(Q("0.02623909"))
    data = json.loads((directory/"ch39-rh-clock-integer-data.json").read_text())["inputs"]
    ref = list(map(cell, data["reference"]))
    p = cell(data["pivot"][0])
    weights = [cell(v)/(Bounds(1)-p) for v in data["rho"]]
    step = 2**100
    @lru_cache(None)
    def value(t):
        terms = [(c*w)/(c+Bounds(t)) for c,w in zip(ref,weights)]
        return Bounds(Q(sum((v.lo*step)//1 for v in terms),step),
                      Q(-sum((-v.hi*step)//1 for v in terms),step))
    @lru_cache(None)
    def balance(total):
        left, right = Q(0), Q(10000)
        assert (value(left)+value(3*left)).lo>total>(value(right)+value(3*right)).hi
        for _ in range(67):
            mid = (left+right)/2
            out = value(mid)+value(3*mid)
            if out.lo>total:
                left = mid
            elif out.hi<total:
                right = mid
            else:
                break
        return left, right
    # Pairwise best maximum error occurs where the two errors have opposite
    # signs and equal magnitudes. This also covers repeated or reversed arcs.
    errors = {}
    for first,second in product(arcs,repeat=2):
        left,right = balance(first+second)
        fitted = Bounds(value(right).lo,value(left).hi)-Bounds(first)
        assert fitted.lo*fitted.hi>0
        errors[first,second] = Bounds(min(abs(fitted.lo),abs(fitted.hi)),
                                      max(abs(fitted.lo),abs(fitted.hi)))
    best = min(errors,key=lambda pair:errors[pair].hi)
    assert all(errors[best].hi<err.lo for pair,err in errors.items() if pair!=best)
    print("Closest arc pair",best,"minimum uniform error",float(errors[best].lo),float(errors[best].hi))
    print("Shortest and marked arcs, minimum error",float(errors[a,b].lo),float(errors[a,b].hi))
    print("All",len(errors),"ordered pairs certified using rational comparisons.")
    assert best==(Q(33,34),Q(31,34)) and Q("0.00315832335797")<errors[best].lo<=errors[best].hi<Q("0.00315832335798")
    assert Q("0.03169364459920")<errors[a,b].lo<=errors[a,b].hi<Q("0.03169364459922")
    robust_pairs(arcs, value, balance, best)



def robust_pairs(arcs, value, balance, best):
    # Mirror the interpolation poles across zero. Every mirrored factor
    # (1-s*z)/(1+s*z) has modulus <= 1 for all nonnegative z.
    from math import prod
    nodes = [j*j for j in range(1,42)]
    factors = {s:prod(Q(s+v,s-v) for v in nodes if s!=v) for s in nodes}
    @lru_cache(None)
    def read(t):
        if t in nodes:
            radius = Q(1,2**96)
        else:
            k = prod((t-s)/(t+s) for s in nodes)
            amp = sum(abs(2*s/(t-s)*k*factors[s]) for s in nodes)
            radius = 2*abs(k)+amp/2**96
        v = value(t)
        return Bounds(v.lo-radius,v.hi+radius)
    enclosures = {}
    for a,b in product(arcs,repeat=2):
        left,right=balance(a+b)
        t=(left+right)/2
        x,y=read(t)-Bounds(a),read(3*t)-Bounds(b)
        # Either side of this one t, monotonicity forces at least the smaller
        # of the two oppositely signed discrepancies. At t both are enclosed.
        if x.hi<0 and y.lo>0:
            x,y=y,x
        assert x.lo>0 and y.hi<0
        enclosures[a,b]=Bounds(min(x.lo,-y.hi),max(x.hi,-y.lo))
    assert all(enclosures[best].hi < box.lo for pair,box in enclosures.items() if pair!=best)
    for pair in (best,(Q(33,34),Q(14,17))):
        v=enclosures[pair]
        print("All positive responses at 96-bit input precision:",pair,float(v.lo),float(v.hi))
    assert enclosures[best].lo>Q("0.00315")


if __name__ == "__main__":
    main()
