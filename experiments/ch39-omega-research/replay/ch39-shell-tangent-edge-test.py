"""Finite simultaneous state/flux probes for the relative-error edge realization bound."""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
from runpy import run_path

from sympy import Matrix

directory=Path(__file__).parent
orient=run_path(str(directory/"ch39-shell-nine-phase-probes.py"))["orientation"]


def main():
    edges=orient()
    D=Matrix.zeros(48,24)
    for e,(a,b) in enumerate(edges):
        D[e,a],D[e,b]=-1,1
    tree=list(D.T.rref()[1])
    chords=[i for i in range(48) if i not in tree]
    right=D[tree,1:].inv()
    bases=[]
    for i in range(48):
        v=[0]+([int(x) for x in right[:,tree.index(i)]] if i in tree else [0]*23)
        twice=[int(e==i)-(v[edges[e][1]]-v[edges[e][0]]) for e in chords]
        bases.append((v,twice))
        check=[v[b]-v[a] for a,b in edges]
        for e,c in zip(chords,twice):
            check[e]+=c
        assert check==[int(e==i) for e in range(48)]
    assert max(abs(x) for v,w in bases for x in v)<=1
    assert max(abs(x) for v,w in bases for x in w)<=2
    # For a sum of up to two different basis vectors, |v_i|<=2, |X_e|<=2;
    # even with |X_e|<=5/2 the Cayley edge residual differs from i*h*r
    # by at most 113*h^2 in Euclidean norm.
    assert 5*2*(Q(25,4)+5)<113
    for bits,h_bits,separating_bits,min_modes in ((104,125,93,14),(179,200,169,21)):
        e=Q(1,2**bits)
        h=Q(1,2**h_bits)
        M=1024
        assert M>=113
        assert (Q(3,2)+113*h)**2<3 and (1+e)*3<4
        reading_error=3*e+(3+113*h)*113*h+(4+M*h)*M*h
        assert reading_error<4*e
        # In each normalized tangent Gram, diagonals are read directly and
        # each offdiagonal is one polarized pair sum divided by two.
        gram_error=(1+47*Q(3,2))*reading_error
        whitened=2*gram_error/(1-gram_error)
        assert whitened<Q(1,2**separating_bits)
        print("Relative input bits",bits,"Cayley exponent",h_bits,"imply at least",48*min_modes,"states.")
    print("Number of nonflat-or-state probes:",48+len(list(combinations(range(48),2))))


if __name__=="__main__":
    main()
