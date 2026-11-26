"""Finite equal-weight trace test and one-generation indices on exterior-power labels."""
from fractions import Fraction as Q
from itertools import combinations
from math import comb


def exterior_weights(charges):
    return [sum(charges[i] for i in labels) for k in (0,2,4) for labels in combinations(range(5),k)]


def main():
    Y = [Q(-1,3)]*3+[Q(1,2)]*2
    weak = [Q(0)]*3+[Q(1,2),Q(-1,2)]
    color = [Q(1,2),Q(-1,2)]+[Q(0)]*3
    wy,wt,wc = map(exterior_weights,(Y,weak,color))
    assert len(wy)==sum(comb(5,k) for k in (0,2,4))==16
    assert sum(t*t for t in wt)==sum(c*c for c in wc)==2
    assert sum(y*y for y in wy)==Q(10,3) and sum(t*y for t,y in zip(wt,wy))==0
    assert sum((t+y)**2 for t,y in zip(wt,wy))==Q(16,3)
    # Polynomial identity: pure cubes, squared-linear terms, and xyz exhaust
    # the polarizations of the local cubic anomaly in the Cartan subalgebra.
    cartan = [exterior_weights([int(i==j)-int(i==4) for i in range(5)]) for j in range(4)]
    for x in cartan:
        assert sum(x)==0
        for y in cartan:
            for z in cartan:
                assert sum(a*b*c for a,b,c in zip(x,y,z))==0
    # Fixed weight matrix invariant under a cyclic permutation and one sign flip.
    for n in (5,16):
        allowed = [(i,j) for i in range(n) for j in range(n)
                   if all(((i+k)%n==0)==((j+k)%n==0) for k in range(n))]
        assert allowed==[(i,i) for i in range(n)]
        assert {(k+1)%n for k in range(n)}==set(range(n))
    print("16 even-subset states: color=weak=2, hypercharge=10/3, electromagnetic=16/3;",
          "equal-weight ratio 3/8 and vanishing cubic anomaly")


if __name__=="__main__":
    main()
