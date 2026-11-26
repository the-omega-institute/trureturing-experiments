"""Exact certificates excluding rank nine and realizing rank fifteen at cutoff 56."""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import fmpq, fmpq_mat

prev = run_path(str(Path(__file__).with_name("ch39-local-memory-rank.py")))
sub = prev["submatrix"]


def verify():
    _, _, _, g0, g1, g2 = prev["MEM"]["coefficients"]()
    delta = g2-g1*g0.solve(g1)
    points = sorted(prev["LOCAL"]["shell_grid"](1)|{(17,17,17)})
    def distance(i, j):
        return prev["LOCAL"]["squared"](points[i], points[j])
    sets = [[4,5,6,7,8,9,10,11,12], [4,6,8,10,11,13,14,15,16], [5,7,9,11,12,13,14,15,16]]
    nulls = []
    for rows in sets:
        cols = sorted(set(range(19))-set(rows))
        assert all(distance(i,j)>56 for i in rows for j in cols)
        fixed = sub(delta, rows, cols)
        num, _ = fixed.numer_denom()
        ker, dim = num.nullspace()
        assert dim == 1 and fixed.rank() == 9
        v = [0]*19
        for j,i in enumerate(cols):
            v[i] = ker[j,0]
        nulls.append(v)
    k = fmpq_mat(nulls)
    assert [i for i in range(19) if distance(0,i)<=56] == [0,2]
    known = k*delta
    witness = fmpq_mat([[k[t,0],k[t,2],known[t,0]] for t in range(3)])
    assert witness.det() != 0

    c = [4,6,8,10]
    r = sorted(set(range(19))-set(c))
    block, cross = sub(delta,r,r), sub(delta,r,c)
    assert all(sub(block,range(j),range(j)).det()>0 for j in range(1,len(r)+1))
    inverse = block.inv()
    x = inverse*cross
    h = sub(delta,c,c)-cross.transpose()*x
    assert [(i,j) for i,j in combinations(range(4),2) if distance(c[i],c[j])>56] == [(0,3)]
    l = r.index(0)
    gain = h[0,3]/(x[l,0]*x[l,3])
    denominator = 1+gain*inverse[l,l]
    assert gain > 0 and denominator > 0
    reduced = fmpq_mat(block)
    reduced[l,l] -= gain/denominator
    corner = cross.transpose()*reduced.solve(cross)
    w = fmpq_mat(delta)
    w[0,0] = reduced[l,l]
    for i,u in enumerate(c):
        for j,v in enumerate(c):
            w[u,v] = corner[i,j]
    assert w.rank() == 15
    assert all(w[i,j]==delta[i,j] for i,j in combinations(range(19),2) if distance(i,j)>56)
    norm = max(sum(abs(v) for v in row) for row in w.tolist())
    inverse_norm = max(sum(abs(v) for v in row) for row in g0.inv().tolist())
    assert norm < fmpq(878,10**9) and inverse_norm < 315000
    assert inverse_norm*norm < (1-fmpq(98,1000))**2
    print("Three forced null vectors contradict row 0 at rank nine.")
    print("Rational completion rank", w.rank(), "row-sum bound", float(norm),
          "auxiliary dimension interval", 29, 34)


if __name__ == "__main__":
    verify()
