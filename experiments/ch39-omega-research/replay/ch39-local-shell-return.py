"""Local directed returns on the integer and denominator-three Behrend spheres."""
from collections import Counter
from math import isqrt
from pathlib import Path
from runpy import run_path

from flint import acb, ctx

GEOMETRY = run_path(str(Path(__file__).with_name("ch39-intrinsic-shell-comparisons.py")))
squared, rotate = GEOMETRY["squared"], GEOMETRY["rotate"]


def shell_grid(den):
    target = 867*den**2
    pts = set()
    for i in range(isqrt(target)+1):
        for j in range(isqrt(target-i*i)+1):
            k = isqrt(target-i*i-j*j)
            if i*i+j*j+k*k == target:
                pts.add((i,j,k))
    pts.remove((17*den,)*3)
    return pts


def scaled(p, factor):
    return tuple(x*factor for x in p)


def graph(den, cutoff):
    nodes = shell_grid(den)
    one, reflected, corner = [scaled(p,den) for p in ((7,17,23),(17,7,23),(11,11,25))]
    exponents = {}
    for _ in range(3):
        for source, charge in ((one,1),(reflected,-1)):
            exponents[source,corner] = charge
            exponents[corner,source] = -charge
        one, reflected, corner = map(rotate,(one,reflected,corner))
    return {p: [(q,exponents.get((p,q),0)) for q in sorted(nodes) if 0<squared(p,q)<=cutoff*den**2]
            for p in nodes}


def routes(adjacency, starts, rotating, steps):
    paths, minima = {}, {}
    for source in sorted(starts):
        target = rotate(source) if rotating else source
        layers = [{(source,0):None}]
        for n in range(1,steps+1):
            next_layer = {}
            for u,k in layers[-1]:
                for v,w in adjacency[u]:
                    next_layer.setdefault((v,k+w),(u,k))
            layers.append(next_layer)
            if (target,1) in next_layer:
                minima.setdefault(source,n)
        if (target,1) in layers[-1]:
            trace, key = [], (target,1)
            for layer in reversed(layers):
                trace.append(key[0])
                key = layer[key]
            path = list(reversed(trace))
            assert path[0]==source and path[-1]==target
            assert sum(dict(adjacency[u])[v] for u,v in zip(path,path[1:])) == 1
            paths[source] = path
    return minima, paths


def verify():
    original = shell_grid(1)
    outer = {p for p in original if sorted(p)==[1,5,29]}
    for p in outer:
        assert min(squared(p,q) for q in original-outer) == 152
        assert sum(squared(p,q)<=151 for q in original if q!=p)==1
    first = (13,13,23)
    near = (first,rotate(first),rotate(rotate(first)))
    def count(p):
        signs = 4 if len(set(p))==2 else 6
        return signs*signs*(3 if p in near else 1)
    small = graph(1,56)
    short, _ = routes(small,original,False,3)
    assert set(short)==original-outer and sum(count(p) for p in short)==408
    assert all(rotate(p) not in original_component(small,p) for p in original)
    assert shell_grid(2)=={scaled(p,2) for p in original}
    for den, cutoff, closed, rotating in ((1,152,5,8),(3,56,6,9)):
        adjacency = graph(den,cutoff)
        sources = {scaled(p,den) for p in original}
        for use_rotation, step in ((False,closed),(True,rotating)):
            minima, paths = routes(adjacency,sources,use_rotation,step)
            print("grid",den,"cutoff",cutoff,"rotate",use_rotation,
                  "minima",dict(Counter(minima.values())),"synchronous",step)
            assert set(paths)==sources and max(minima.values())==step
            if use_rotation and den==3:
                path = paths[(3,15,87)]
                print("Nine-step witness, numerators over three",path)
    assert len(shell_grid(3))==75
    ctx.prec = 100
    for n in range(11485,11491):
        rho = acb.zeta_zero(n)
        w = (rho-1)/rho
        assert all(not (w**m-1).contains(0) for m in range(1,11))


def original_component(adjacency, start):
    component, pending = {start}, [start]
    while pending:
        for v,_ in adjacency[pending.pop()]:
            if v not in component:
                component.add(v)
                pending.append(v)
    return component


if __name__ == "__main__":
    verify()
