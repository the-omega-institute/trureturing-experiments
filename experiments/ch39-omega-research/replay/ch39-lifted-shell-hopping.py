"""Exact intertwiners on the nineteen-point shell and its faithful nine-cycle lift."""

from collections import Counter
from pathlib import Path
from runpy import run_path

import sympy as sp


def verify():
    local = run_path(str(Path(__file__).with_name("ch39-local-product-contact.py")))
    sites, _ = local["lifted_sites"]()
    graph, _ = local["lifted_laplacian"]()
    old = sorted(set(sites))
    sizes = Counter(sites)
    indices = {p: old.index(p) for p in old}
    owners = [indices[p] for p in sites]
    n, k = len(old), len(sites)
    graph = sp.Matrix(graph.tolist())
    adjacency = sp.diag(*graph.diagonal()) - graph
    projection = sp.zeros(k,n)
    for i, x in enumerate(owners):
        projection[i,x] = 1
    multiplicities = projection.T * projection
    normalized = projection * sp.diag(*(1/sp.sqrt(m) for m in multiplicities.diagonal()))
    assert normalized.T * normalized == sp.eye(n)

    base = sp.zeros(n,n)
    distance = local["geometry"]["squared"]
    for i, x in enumerate(old):
        for j, y in enumerate(old):
            base[i,j] = int(0 < distance(x,y) <= 152)
    assert projection * base * projection.T == adjacency
    kinds = Counter(sizes[old[i]] * sizes[old[j]] for i in range(n)
                    for j in range(i+1,n) if base[i,j])
    assert kinds == {1:27, 3:12}
    assert sum(adjacency) == 2*63
    degree = list(base*sp.ones(n,1))
    triple = [i for i,p in enumerate(old) if sizes[p]==3]
    assert len(triple)==3 and all(degree[i]==4 for i in triple)

    rotate = lambda p: (p[1],p[2],p[0])
    first = (13,13,23)
    third = rotate(rotate(first))
    lift = sp.zeros(k)
    for i,p in enumerate(sites):
        tag = i-sites.index(p)
        target = sites.index(rotate(p)) + (tag+int(p==third)) % sizes[p]
        lift[target,i] = 1
    rotation = sp.zeros(n)
    for i,p in enumerate(old):
        rotation[indices[rotate(p)],i] = 1
    assert lift*projection == projection*rotation
    assert lift**9 == sp.eye(k) and lift**3 != sp.eye(k)
    assert rotation*base == base*rotation

    rates = adjacency * sp.diag(*(sp.Rational(1,sizes[p]) for p in sites))
    generator = rates - sp.diag(*(rates*sp.ones(k,1)))
    base_laplacian = sp.diag(*degree)-base
    mass = sp.diag(*(sp.Rational(1,n*sizes[p]) for p in sites))
    assert mass*generator == generator.T*mass
    assert generator*projection == -projection*base_laplacian
    assert lift*generator == generator*lift
    contrasts = []
    for i in triple:
        fiber = [k for k,x in enumerate(owners) if x==i]
        for j in fiber[1:]:
            col = sp.zeros(len(sites),1)
            col[fiber[0]], col[j] = 1,-1
            contrasts.append(col)
    dark = sp.Matrix.hstack(*contrasts)
    assert dark.rank()==6 and projection.T*dark==sp.zeros(n,6)
    assert generator*dark == -4*dark

    hopping = normalized*base*normalized.T
    laplacian = sp.diag(*(degree[x] for x in owners))-hopping
    assert hopping*normalized == normalized*base
    assert laplacian*normalized == normalized*base_laplacian
    assert lift*hopping == hopping*lift
    assert hopping*dark == sp.zeros(k,6) and laplacian*dark == 4*dark
    variable = sp.symbols("t")
    assert sp.cancel(lift.charpoly(variable).as_expr()
                     /rotation.charpoly(variable).as_expr()) == variable**6+variable**3+1
    visible = normalized.T*adjacency*normalized
    assert {visible[i,j] for i in range(n) for j in range(i+1,n) if base[i,j]} == {1,sp.sqrt(3)}

    classical_scale = sp.Rational(1,2)
    assert (1-classical_scale) == 3*classical_scale-1
    amplitude = 2/(1+sp.sqrt(3))
    assert sp.simplify(1-amplitude-(2-sp.sqrt(3)))==0
    assert sp.simplify(sp.sqrt(3)*amplitude-1-(2-sp.sqrt(3)))==0
    print("Base edges: 27 singleton-singleton, 12 singleton-triple; lifted edges: 63.")
    print("Six contrasts: generator -4, adjacency 0, normalized Hamiltonian lift of Laplacian 4.")
    print("Equal coefficient errors in individual rates / visible matrix entries: 1/2, 2-sqrt(3).")


if __name__ == "__main__":
    verify()
