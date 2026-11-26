"""Spatial range at which a ninth clock can hide a new visible relay mode."""
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from runpy import run_path

from flint import fmpq_mat

prev = run_path(str(Path(__file__).with_name("ch39-local-memory-rank.py")))
sub = prev["submatrix"]


def verify():
    d, pair, *_ = prev["MEM"]["coefficients"]()
    b = fmpq_mat([row[:19] for row in pair.tolist()])
    square = d**2*b
    space = fmpq_mat([x+y for x,y in zip(pair.tolist(), square.tolist(), strict=True)])
    num, _ = space.transpose().numer_denom()
    kernel, size = num.nullspace()
    assert size == 4
    normal = sub(fmpq_mat(kernel), range(57), range(4))
    n, rank = (normal.transpose()*d*square).rref()
    assert rank == 4
    points = sorted(prev["LOCAL"]["shell_grid"](1)|{(17,17,17)})
    dist = prev["LOCAL"]["squared"]
    left, right = [0,1,2,4], [18,17,15,12]
    assert min(dist(points[i],points[j]) for i in left for j in right) == 296
    assert all(sub(n,range(4),part).det() for part in (left, right))
    symmetric = list(combinations_with_replacement(range(4), 2))
    constraints = []
    for i,j in combinations(range(19), 2):
        if dist(points[i],points[j]) >= 512:
            constraints.append([n[a,i]*n[b,j] + (a!=b)*n[b,i]*n[a,j] for a,b in symmetric])
    assert fmpq_mat(constraints).rank() == 10
    witness = [n[3,i] for i in range(19)]
    support = [points[i] for i,x in enumerate(witness) if x]
    assert len(support) == 6 and {tuple(sorted(p)) for p in support} == {(7,17,23)}
    assert sorted(x for x in witness if x) == [-1,-1,-1,1,1,1]
    assert max(dist(p,q) for p,q in combinations(support,2)) == 512
    print("Nine clocks uniquely determine all frequencies exactly for cutoffs below 512.")
    print("Seven clocks require all 57 relay states for cutoffs below 296.")


if __name__ == "__main__":
    verify()
