"""Exact spatial range of a conserved family charge on the 25-site shell."""
from collections import defaultdict
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import fmpq_mat as Matrix

leads = run_path(str(Path(__file__).with_name("ch39-unit-lead-gauge-shell.py")))
local = leads["local"]
squared = local["geometry"]["squared"]


def verify():
    points, pole = local["lifted_sites"]()
    sites = [(p, points[:i].count(p)) for i, p in enumerate(points)]
    edges, _ = local["lifted_laplacian"]()
    laplacian = Matrix(edges)
    n = len(sites)
    probes = sorted({p for p in points if squared(p, pole) == 68})
    charge = Matrix([[int(i == j)-1 for j in range(3)] for i in range(3)])
    symmetries = [[sites.index((tuple(p[a] for a in perm), k)) for p, k in sites]
                  for perm in permutations(range(3))]
    rotation = [sites.index(((p[1], p[2], p[0]), k)) for p, k in sites]
    shift = Matrix([[int(rotation[j] == i) for j in range(n)] for i in range(n)])
    conserved = -shift-shift.transpose()
    identity = Matrix([[int(i == j) for j in range(n)] for i in range(n)])
    assert conserved*laplacian == laplacian*conserved
    assert (conserved-identity)*(conserved+2*identity) == 0*identity

    hops = [[0 if i == j else (1 if edges[i,j] else n) for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                hops[i][j] = min(hops[i][j], hops[i][k]+hops[k][j])
    pairs = [(i,j) for i in range(n) for j in range(i,n)]
    support = [(i,j) for i,j in pairs if conserved[i,j]]
    assert max(hops[i][j] for i,j in support) == 4
    assert max(squared(points[i],points[j]) for i,j in support) == 1376

    orbits = defaultdict(list)
    for i,j in pairs:
        label = min(tuple(sorted((p[i],p[j]))) for p in symmetries)
        orbits[label].append((i,j))
    for copies in (1,2,3):
        common, contrast = leads["terminal_modes"](copies)
        assert contrast > common > 0
        response = Matrix([[(common-contrast)/3+int(i==j)*contrast for j in range(3)]
                           for i in range(3)])
        eye = Matrix([[int(i==j) for j in range(3)] for i in range(3)])
        assert 3*response-(common+2*contrast)*eye == (contrast-common)*charge
        ports = Matrix([[int(p==q and k<copies) for q in probes] for p,k in sites])
        assert conserved*ports == ports*charge
        commutator_zeros = [0]*(n*(n-1)//2)
        target = commutator_zeros + (ports*charge).entries()
        chord_rows, hop_rows = [], []
        for (i,j), orbit in orbits.items():
            basis = Matrix(n,n)
            for a,b in orbit:
                basis[a,b] = basis[b,a] = 1
            bracket = basis*laplacian-laplacian*basis
            coeff = [bracket[a,b] for a,b in pairs if a<b] + (basis*ports).entries()
            if squared(points[i],points[j]) < 1376:
                chord_rows.append(coeff)
            if hops[i][j] < 4:
                hop_rows.append(coeff)
        ranks = []
        for rows in (chord_rows,hop_rows):
            rank = Matrix(rows).rank()
            assert Matrix(rows+[target]).rank() == rank+1
            ranks.append(rank)
        assert ranks == ([57,56] if copies==3 else [61,60])
        print(copies, "copies: rank obstruction", ranks, "sharp squared chord 1376, graph distance 4")


if __name__ == "__main__":
    verify()
