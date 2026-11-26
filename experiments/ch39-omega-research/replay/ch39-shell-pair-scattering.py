"""Finite shell registers and the exact signed multiplicity of a symmetric reflection."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import arb,acb,ctx

ROOT = Path(__file__).parent
UNIT = run_path(str(ROOT/"ch39-unitary-shell-phase.py"))
PHASE,a = UNIT["PHASE"],UNIT["a"]


def channels():
    shell = [v for v in product(range(30),repeat=3) if sum(x*x for x in v)==867]
    orbits = Counter(tuple(sorted(v)) for v in shell)
    assert sorted(orbits.values())==[1,3,3,6,6]
    pole, = [v for v in shell if len(set(v))==1]
    def dist(x,y):
        return sum((i-j)**2 for i,j in zip(x,y))
    pairs = [(x,y) for x,y in combinations(shell,2) if sorted(x)!=sorted(y)]
    near = min(dist(pole,y) for y in shell if y!=pole)
    kept = [(x,y) for x,y in pairs if not (pole in (x,y) and dist(x,y)==near)]
    assert len(pairs)==135 and near==68 and len(kept)==132
    assert all(sorted(y if x==pole else x)==[13,13,23] for x,y in set(pairs)-set(kept))
    return len(kept)


def locked_gains(d,h,gamma):
    def stretch(t):
        u = d/t
        return (1+u*u+u*(2+u*u).sqrt()).sqrt()
    outputs = []
    for sign in (1,-1):
        lo,hi = F(".6"),F("1.2")
        for _ in range(32):
            mid = (lo+hi)/2
            t = a(mid)
            gap = h*t**4*stretch(t)**sign-d
            assert gap>0 or gap<0
            if gap>0: hi=mid
            else: lo=mid
        t = a(lo).union(a(hi))
        gain = gamma/((t/d)**4-a("1/4")).sqrt()
        lower_integer = int(gain.floor().unique_fmpz())
        assert lower_integer < gain < lower_integer+1
        outputs.append(gain)
    return outputs


def verify():
    ctx.prec = 150
    nslots = channels()
    for n,kind,r,v,state in PHASE["MATCH"]["matched_points"]():
        zero = acb.zeta_zero(n)
        w = (zero-1)/zero
        z = 1/state[0]
        for gain in (nslots,19**2):
            # With one external input of 1, the vector of internal incoming waves is constant.
            internal = 2*w/(gain+1-(gain-1)*w)
            average = (1+gain*internal)/(gain+1)
            assert (w*(2*average-internal)-internal).contains(0)
            reflected = 2*average-1
            assert (abs(reflected)-1).contains(0)
            UNIT["certify_masses"](reflected,r,z)
        _,minimum,_ = PHASE["phase_bounds"](r,z)
        assert 1/((zero.imag/19)**2+a("1/4")).sqrt()<minimum
        # Every feasible effective height is greater than 15.
        assert zero.imag/19**3 < 15
        d,h = [(a(x)*z).sqrt()*r for x in ("13/10","7/2")]
        first,second = locked_gains(d,h,zero.imag)
        print(n,kind,"aligned gains",first.str(10),second.str(10))


if __name__=="__main__":
    verify()
