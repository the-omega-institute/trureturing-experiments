"""One- and two-endpoint observation ranks and a finite modular reading."""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from sympy import Matrix, zeros

model = run_path(str(Path(__file__).with_name("ch39-finite-relative-yukawa.py")))


def main():
    pairs = list(combinations(range(4),2))
    rows = [(p,c) for p in pairs for c in (-1,1)]
    points = list(permutations((0,1,3,4)))
    a,b = points[:2]
    dot = lambda p,s:sum(x*y for x,y in zip(p,s))
    records = [(dot(p,a),dot(p,b)) for p in points]
    _,up,down = model["benchmark"]()
    target = Matrix([-up[p]-c*down[p].conjugate() for p,c in rows])
    for width, expected in [(1,(13,9)),(2,(24,12))]:
        keys = sorted({r[:width] for r in records})
        counts = zeros(12,len(keys))
        for p,r in zip(points,records):
            root = tuple(sorted((p.index(3),p.index(4))))
            chi = (-1)**sum(p[i]>p[j] for i,j in pairs)
            counts[rows.index((root,chi)),keys.index(r[:width])]+=1
        assert (len(keys),counts.rank())==expected
        if width == 1:
            assert counts.row_join(target).rank()==10
        else:
            assert len(keys)==len(set(records))==24
    def valid(m,i,j):
        return len({(i*s+j*t)%m for s,t in records})==24
    # All smaller moduli from 24 upward; below 24 the pigeonhole principle applies.
    assert all(not valid(m,i,j) for m in range(24,39) for i in range(m) for j in range(m))
    assert valid(39,1,7)
    print("One-anchor rank 9 with incompatible benchmark; two-anchor rank 12, all 24 vertices distinguished.")
    print("The integer residue (a.p + 7 b.p) mod 39 distinguishes the whole shell.")


if __name__=="__main__":
    main()
