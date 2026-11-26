"""Finite switch geometry, rational reflection tests, and exact fluctuation checks."""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import comb
from pathlib import Path
from runpy import run_path

import sympy as sp

geometry = run_path(str(Path(__file__).with_name("ch39-regge-volume-reflection.py")))
vertex, lift, reflect = (geometry[name] for name in ("vertex","lift","reflect"))


def switch_geometry():
    pos = set()
    for start in product(range(4),repeat=4):
        if start[0] in (0,1):
            x = list(lift(start))
            x[0] += 1
            pos.add(frozenset((start,vertex(x))))
    neg = {frozenset(map(reflect,edge)) for edge in pos}
    assert len(pos) == len(neg) == 128 and pos.isdisjoint(neg)
    degree = Counter()
    for start in product(range(4),repeat=4):
        for order in permutations(range(4)):
            x = list(lift(start))
            path = [start]
            for i in order:
                x[i] += 1
                path.append(vertex(x))
            selected = [frozenset((u,v)) for i,u in enumerate(path) for v in path[i+1:]
                        if frozenset((u,v)) in pos or frozenset((u,v)) in neg]
            assert len(selected) == 1
            degree[selected[0]] += 1
    assert degree == {edge:24 for edge in pos|neg}
    return len(pos|neg)


def moments(count,weight):
    masses = [comb(count,k)*weight(k) for k in range(count+1)]
    total = sum(masses)
    mean = Fraction(sum(k*m for k,m in enumerate(masses)),total)
    variance = sum((k-mean)**2*m for k,m in enumerate(masses))/total
    parity = Fraction(sum((-1)**k*m for k,m in enumerate(masses)),total)
    return mean,variance,parity


def verify():
    count = switch_geometry()
    for bias in (Fraction(1,3),Fraction(1),Fraction(5,2)):
        mean,var,_ = moments(count,lambda k:bias**k)
        p = bias/(1+bias)
        assert mean == count*p and var == count*p*(1-p)
    mean,var,parity = moments(count,lambda k:1+Fraction((-1)**k,3))
    assert (mean,var,parity) == (128,64,Fraction(1,3))
    assert Fraction(2,3840) == Fraction(1,1920)
    r = sp.symbols("r")
    kink = sp.Matrix([[1,1,1],[1,1,r],[1,r,r**2]])
    assert sp.factor(kink.det()+(r-1)**2) == 0
    print("Switches 256, incident simplices 24 each; binomial variance and parity check:", mean,var,parity)
    print("Three half-volumes at a log-slope change:", sp.factor(kink.det()))


if __name__ == "__main__":
    verify()
