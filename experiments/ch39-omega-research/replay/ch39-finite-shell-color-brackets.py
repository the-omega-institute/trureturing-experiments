"""Finite rational bracket tables of the shortest-arc records in the doubled tangent space."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

cut = run_path(str(Path(__file__).with_name("ch39-finite-permutation-arc-cut.py")))
mul,eye,rank,walk = (cut[name] for name in ("mul","eye","rank","walk"))


def add(a,b,s=1):
    return [[x+s*y for x,y in zip(row,col)] for row,col in zip(a,b)]


def bracket(a,b):
    return add(mul(a,b),mul(b,a),-1)


def block(a,b):
    return [u+[-x for x in v] for u,v in zip(a,b)]+[v+u for u,v in zip(a,b)]


def flat(m):
    return [x for row in m for x in row]


def generated(generators):
    basis = list(generators)
    sizes = [len(basis)]
    for _ in range(6):
        snapshot = basis.copy()
        for b in snapshot:
            for g in generators:
                term = bracket(g,b)
                if rank([flat(t) for t in basis]+[flat(term)])>len(basis):
                    basis.append(term)
        if len(basis)==len(snapshot):
            return basis,sizes
        sizes.append(len(basis))
    raise AssertionError("unexpected bracket depth")


def main():
    square = [(0,1),(3,4)]*2
    hexagon = [(0,1),(1,3)]*3
    a,ai,b,bi = (walk(s) for s in (square,square[::-1],hexagon,hexagon[::-1]))
    assert mul(a,ai)==mul(b,bi)==eye(3)
    skew = [add(a,ai,-1),add(b,bi,-1)]
    assert len(generated(skew)[0])==3
    even = add(a,ai)
    trace = sum(even[i][i] for i in range(3))
    for i in range(3):
        even[i][i] -= Q(trace,3)
    zero = [[0]*3 for _ in range(3)]
    generators = [block(h,zero) for h in skew]+[block(zero,even)]
    basis,sizes = generated(generators)
    metric = [[1,0,0],[0,10,12],[0,12,17]]
    metric6 = block(metric,zero)
    for H in basis:
        assert add(mul(list(zip(*H)),metric6),mul(metric6,H)) == [[0]*6 for _ in range(6)]
        assert sum(H[i][i] for i in range(3))==sum(H[i+3][i] for i in range(3))==0
    assert len(basis)==8
    # Commuting real tangent endomorphisms form a one-dimensional subspace.
    units = [[[int(i==r and j==c) for j in range(3)] for i in range(3)]
             for r in range(3) for c in range(3)]
    assert rank([flat(bracket(m,skew[0]))+flat(bracket(m,skew[1])) for m in units])==8
    print("Real rotation bracket rank 3; doubled bracket ranks",sizes,
          "; scalar centralizer on the original tangent triple.")


if __name__=="__main__":
    main()
