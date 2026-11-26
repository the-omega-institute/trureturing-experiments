"""How the spatial range of an old-shell counterterm limits three-clock compression."""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import fmpq, fmpq_mat

MEM = run_path(str(Path(__file__).with_name("ch39-three-clock-minimal-memory.py")))
LOCAL = MEM["REFINE"]["LOCAL"]


def submatrix(matrix, rows, columns):
    return fmpq_mat([[matrix[i,j] for j in columns] for i in rows])


def verify():
    _, _, _, g0, g1, g2 = MEM["coefficients"]()
    shift = g2-g1*g0.solve(g1)
    points = sorted(LOCAL["shell_grid"](1)|{(17,17,17)})
    pairs = {(i,j):LOCAL["squared"](p,q) for i,p in enumerate(points) for j,q in enumerate(points)}
    assert all(shift[i,j]!=0 for i,j in combinations(range(19),2))
    assert max(pairs.values())==1568
    far = [(i,j) for i,j in combinations(range(19),2) if pairs[i,j]==1568]
    b, = {shift[i,j] for i,j in far}
    a, = {shift[i,j] for i,j in combinations(range(19),2) if pairs[i,j]==1376}
    assert b<0 and a<0 and abs(b)<-2*a
    vertices = {i for ij in far for i in ij}
    assert len(far)==3 and len(vertices)==6
    # Two triangles of 120-degree-separated vectors; far pairs form the matching.
    seed = min(vertices)
    triple = {seed}|{j for j in vertices if pairs[seed,j]==1376}
    other = vertices-triple
    for sides in (triple,other):
        assert len(sides)==3 and all(pairs[i,j]==1376 for i,j in combinations(sides,2))
    assert all((i in triple)!=(j in triple) for i,j in far)
    assert all(d in (1376,1568) for d in pairs.values() if d>1152)
    assert max(sum(abs(x) for x in row) for row in g0.inv().tolist())<315000
    assert 315000*6*(-2*a)<fmpq(1,5000)
    # Separated sets of old components give a forced rank-nine rectangular block.
    part = [4,5,6,7,8,9,10,11,12]
    rest = sorted(set(range(19))-set(part))
    assert all(pairs[i,j]>56 for i in part for j in rest)
    assert submatrix(shift,part,rest).rank()==9
    clique = next(c for c in combinations(range(19),3) if all(pairs[i,j]<=56 for i,j in combinations(c,2)))
    rest = sorted(set(range(19))-set(clique))
    cross = submatrix(shift,rest,clique)
    correction = submatrix(shift,clique,clique)-cross.transpose()*submatrix(shift,rest,rest).solve(cross)
    residual = fmpq_mat(shift)
    for i, u in enumerate(clique):
        for j, v in enumerate(clique):
            residual[u,v] -= correction[i,j]
    assert residual.rank()==16
    print("Far entries", float(a),float(b),"exact dimensions by cutoff 1152,1376,1568:",
          21,20,19,"cutoff 56 bounds",28,35)


if __name__ == "__main__":
    verify()
