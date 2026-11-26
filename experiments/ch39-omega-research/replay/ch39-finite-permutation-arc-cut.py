"""Exact first noncommuting arc cutoff on the 24 permutations of (0,1,3,4)."""
from itertools import combinations, permutations, product
from pathlib import Path
from runpy import run_path

folder = Path(__file__).resolve().parent
arc = run_path(str(folder/"ch39-finite-shortest-arc-holonomy.py"))
record = run_path(str(folder/"ch39-finite-behrend-zero-record.py"))
mul,eye,rank = arc["mul"],arc["eye"],arc["rank"]
start = (0,1,3,4)


def walk(swaps):
    p = start
    W = eye(4)
    for i,j in swaps:
        dest = tuple(j if x==i else i if x==j else x for x in p)
        W = mul(arc["link"](p,dest),W)
        p = dest
    assert p==start
    F = [[1,0,0],[0,-3,-4],[0,1,0],[0,0,1]]
    B = mul(W,F)
    V = [B[0],B[2],B[3]]
    assert mul(F,V)==B
    return V


def main():
    pts = list(permutations(start))
    cols,weights = record["basis_and_weights"]()
    # Represent four-coordinate shell vectors in the orthonormal record basis.
    def embed(p,c):
        return tuple(sum(v*z[k] for v,z in zip(p,c)) for k in (0,1))
    for p,q in product(pts,repeat=2):
        real = sum(w*record["dot"](embed(p,c),embed(q,c)) for w,c in zip(weights,cols))
        assert real==sum(x*y for x,y in zip(p,q))
    edges = [(p,q,26-sum(x*y for x,y in zip(p,q))) for p,q in combinations(pts,2)]
    for cut,n in [(1,24),(2,36),(3,36),(4,48)]:
        E = [(p,q) for p,q,c in edges if c<=cut]
        assert len(E)==n
        if cut < 4:
            assert all(all((x<=1)==(y<=1) for x,y in zip(p,q)) for p,q in E)
    A = walk([(0,1),(3,4)]*2)
    B = walk([(0,1),(1,3)]*3)
    assert mul(A,B)!=mul(B,A)
    words = layer = [eye(3)]
    for d in (3,7,9):
        layer = [mul(w,g) for w,g in product(layer,(A,B))]
        words = words+layer
        assert rank([[x for row in w for x in row] for w in words])==d
    print("576 exact pair readings. First noncommuting squared-chord/2 cutoff is 4.",
          "Word spans through length 1,2,3 have dimensions 3,7,9.")


if __name__=="__main__":
    main()
