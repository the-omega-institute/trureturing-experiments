"""Exact tests for the Reed-Solomon dual description of all additional mirror pairs."""
from collections import Counter
from fractions import Fraction as F
from itertools import product,combinations
from math import prod
from flint import fmpq_mat

# Distinct rational points make the dual-code test exact without a number-field dependency.
# The dimension argument in the note applies to any six distinct complex points.
z=(0,1,2,4,6,10)
weights=[F(1,prod(v-w for j,w in enumerate(z) if i!=j)) for i,v in enumerate(z)]
assert all(sum(w*v**k for w,v in zip(weights,z))==0 for k in range(5))


def mirror(directions):
    evals=fmpq_mat([[v**k*a for k in range(3) for a in row] for v,row in zip(z,directions)])
    annihilator=[]
    for v,row in zip(z,directions):
        for a,b in combinations(range(3),2):
            annihilator.append([v**k*(row[a]*int(i==b)-row[b]*int(i==a)) for k in range(3) for i in range(3)])
    opposite=6-evals.rank()
    assert 9-fmpq_mat(annihilator).rank()==opposite
    return opposite


def main():
    for degree in range(3):
        for exceptions in range(4-degree):
            dirs=[tuple(v**i if i<=degree else 0 for i in range(3)) for v in z]
            for k in range(exceptions):
                dirs[k]=(1,7+k,17+2*k)
            assert mirror(dirs)==3-degree-exceptions
    basis=[(1,0,0),(0,1,0),(0,0,1)]
    strata=Counter()
    balanced=0
    for ids in product(range(3),repeat=6):
        q=mirror([basis[j] for j in ids])
        assert q==sum(max(0,ids.count(j)-3) for j in range(3))
        strata[q]+=1
        balanced+=all(ids.count(j)==2 for j in range(3))
    assert balanced==90
    assert strata=={0:510,1:180,2:36,3:3}
    print("Cartan assignments by number of mirror pairs",sorted(strata.items()),"balanced",balanced)
    # A primitive degree d map has 3(d+1)-1 parameters, plus two for each omitted site.
    for q in (1,2,3):
        dimensions=[3*(d+1)-1+2*k for d in range(3) for k in range(3) if d+k<=3-q]
        print("at least",q,"mirror pairs: complex codimension",12-max(dimensions))


if __name__=="__main__":
    main()
