"""Test the zero-resolvent threshold against all modes of a Behrend shell arc graph."""
from itertools import combinations, product
from math import isqrt
from pathlib import Path
from runpy import run_path
from flint import arb, fmpz_mat

arith = run_path(str(Path(__file__).with_name("ch39-shell-rh-parent-threshold.py")))
Z = arith["zero_resolvent"]
target = arb(12)/13


def spectrum(radius_squared):
    points = [p for p in product(range(isqrt(radius_squared)+1),repeat=3)
              if sum(x*x for x in p) == radius_squared]
    lengths = {}
    for i,j in combinations(range(len(points)),2):
        distance = radius_squared-sum(a*b for a,b in zip(points[i],points[j]))
        lengths.setdefault(distance,[]).append((i,j))
    lap = fmpz_mat(len(points),len(points))
    components = []
    for distance, edges in sorted(lengths.items()):
        for i,j in edges:
            lap[i,i] += 1
            lap[j,j] += 1
            lap[i,j] -= 1
            lap[j,i] -= 1
        nullity = len(points)-lap.rank()
        components.append((distance,nullity))
        if nullity == 1:
            break
    roots = lap.charpoly().complex_roots()
    assert all(z.imag == 0 and (z.real == 0 or z.real>0) for z,_ in roots)
    return components, [(z.real,n) for z,n in roots]


def increasing_root(f,upper,steps=58):
    lo,hi = arb(0),arb(upper)
    assert f(lo)<0 and f(hi)>0
    for _ in range(steps):
        m = (lo+hi)/2
        v = f(m)
        if v<0:
            lo = m
        else:
            assert v>0
            hi = m
    assert f(lo)<0 and f(hi)>0
    return lo.union(hi)


def ratio_squared(s,mass_squared):
    return Z((s*s+mass_squared).sqrt())/Z(s)


if __name__ == "__main__":
    for r in (17,26,29,867):
        components,roots = spectrum(r)
        q0 = ratio_squared(arb(0),r)
        assert (q0>target) if r != 867 else (q0<target and 2*q0>target)
        print("shell",r,"initial single-mode weight",q0.str(13),"component counts",components)
    assert components == [(6,16),(16,13),(26,7),(28,7),(34,4),(36,4),(68,4),(76,1)]
    def total(s,k):
        base = Z(s)
        return sum(n*Z((s*s+867+k*lam).sqrt()) for lam,n in roots)/base
    s_limit = increasing_root(lambda s: ratio_squared(s,867)-target, 60)
    k_crit = increasing_root(lambda k: target-total(arb(0),k), 100000000)
    print("critical angular radius in degrees",(1-arb(76)/867).acos()*180/arb.pi())
    print("limiting Q/mu at infinite stiffness",s_limit.str(14))
    print("critical stiffness",k_crit.str(14))
    for multiple in (2,10):
        k = multiple*k_crit
        lower = increasing_root(lambda s: total(s,k.lower())-target,60)
        upper = increasing_root(lambda s: total(s,k.upper())-target,60)
        print("Q/mu at stiffness multiple",multiple,lower.union(upper).str(13))
