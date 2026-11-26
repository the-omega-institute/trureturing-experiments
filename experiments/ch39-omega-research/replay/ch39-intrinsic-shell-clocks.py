"""Construct the intrinsic shell clock and its four integer valuation coordinates."""
from itertools import combinations, permutations, product
from math import lcm
from pathlib import Path
from runpy import run_path

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form


def signature(c, d):
    primes = (3, 7, 11, 31)
    roots = (1, 2, 4, 13)
    den = lcm(int(c.q), int(d.q))
    num, imag = int(c*den), int(d*den)
    assert num*num + 17*imag*imag == den*den
    remaining = den
    out = []
    for p, r in zip(primes,roots):
        assert (r*r+17) % p == 0
        count = 0
        while remaining%p == 0:
            remaining //= p
            count += 1
        if count:
            assert ((num+imag*r)%p == 0) != ((num-imag*r)%p == 0)
        out.append(count if (num+imag*r)%p == 0 else -count)
    assert remaining == 1
    return out


def main():
    book = run_path(str(Path(__file__).with_name("ch39-omega-projector-clock.py")))
    points = [sp.Matrix(p) for p in product(range(5),repeat=3) if sum(t*t for t in p)==17]
    eye = sp.eye(4)
    j = book["multiply_by_vector"](points[0], eye, right=True)
    def involution(p):
        return book["multiply_by_vector"](p,eye,right=False)*j/17
    inner = [sp.Matrix(p) for p in ((3,2,2),(2,3,2),(2,2,3))]
    vol = sp.zeros(4)
    for i,k,l in permutations(range(3)):
        inversions = sum(a>b for a,b in combinations((i,k,l),2))
        vol += (-1)**inversions * involution(inner[i])*involution(inner[k])*involution(inner[l])
    assert 289*vol == 42*j
    vectors,types=[],set()
    for p,q,r in combinations(points,3):
        a=(17+p.dot(q)+r.dot(p)+q.dot(r))/68
        b=sp.Matrix.hstack(p,q,r).det()/1156
        c=(a*a-17*b*b)/(a*a+17*b*b)
        d=2*a*b/(a*a+17*b*b)
        vectors.append(signature(c,d))
        types.add((c,abs(d)))
    assert len(types)==13
    lattice=hermite_normal_form(sp.Matrix(vectors).T)
    assert lattice==sp.Matrix([[2,1,1,1],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    print("The antisymmetrized product recovers J. All 84 triangle turns generate a rank-four valuation lattice of index two.")


if __name__ == "__main__":
    main()

