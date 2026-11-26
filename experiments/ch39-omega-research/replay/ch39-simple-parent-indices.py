"""Finite index and dimension checks for a common Pati-Salam gauge parent."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent


def colour(a,b,c):
    rows = (a+b+c,b+c,c,0)
    boxes = sum(rows)
    size = prod(F(rows[i]-rows[j]+j-i,j-i) for i in range(4) for j in range(i+1,4))
    casimir = (4*boxes + sum(s*(s+1-2*i) for i,s in enumerate(rows,1)) - F(boxes**2,4))/2
    index = 2*size*casimir/15
    assert size.denominator == index.denominator == 1
    return int(size),int(index)


def block(labels):
    a,b,c,l,r = labels
    dim,idx = colour(a,b,c)
    left,right = l+1,r+1
    return (dim*left*right, (idx*left*right,
            dim*right*(left**3-left)//6, dim*left*(right**3-right)//6))


def verify():
    source = run_path(str(root/"ch39-rational-gauge-flow.py"))
    y,l,c = source["RATIOS"]
    assert (15*c,15*l,15*(5*y-2*c)/3) == (15,16,14)
    central = F(15*15,15+4)+F(3*16,16+2)+F(3*14,14+2)
    assert central == F(7813,456) and 17 < central < 18
    assert F(9*19,11) < central  # Sp(18) is already too small.
    # Stop at dimension 56, enough to certify the larger target's minimum.
    assert colour(6,0,0)[0] > 56 and colour(0,4,0)[0] > 56
    assert (6**3-6)//6 > 32
    entries = []
    primitive = (15,16,14)
    double = tuple(2*x for x in primitive)
    for a,b,c,l,r in product(range(6),range(4),range(6),range(5),range(5)):
        if c>a:
            continue  # conjugating a colour block preserves both dimension and indices
        dim,idx = block((a,b,c,l,r))
        if dim <= 56 and any(idx) and all(x<=cap for x,cap in zip(idx,double)):
            entries.append((dim,idx,l%2))

    @lru_cache(None)
    def minimum(target, needs_doublet):
        if not any(target):
            return 1000 if needs_doublet else 0
        best = 1000
        for dim,idx,odd in entries:
            rest = tuple(x-y for x,y in zip(target,idx))
            if min(rest) >= 0:
                best = min(best,dim+minimum(rest,needs_doublet and not odd))
        return best

    assert minimum(primitive,False) == 35
    assert minimum(primitive,True) == 41
    assert minimum(double,False) == 56
    witnesses = (
        (35,[(0,0,0,0,2),(1,0,0,2,0),(2,0,0,0,1)]),
        (41,[(0,0,0,0,1)]*2 + [(0,0,0,2,2),(1,0,0,1,0),(1,1,0,0,0)]))
    for expected,labels in witnesses:
        blocks = [block(w) for w in labels]
        assert sum(d for d,_ in blocks) == expected
        assert tuple(map(sum,zip(*(idx for _,idx in blocks)))) == primitive
    # The 41-dimensional representation splits as 20 + 21, with the
    # highest-weight (1,1,0) colour block alone in the first factor.
    split = witnesses[1][1]
    first = block(split[-1])
    second = [block(w) for w in split[:-1]]
    assert first == (20,(13,0,0))
    assert sum(d for d,_ in second) == 21
    assert tuple(map(sum,zip(*(idx for _,idx in second)))) == (2,16,14)
    # For Sp at k=1 all self-conjugate SU(4) constituents with available index
    # are real: 1, the real six and the adjoint. For k=3 the dimension is <=24.
    small = [((a,b,c),colour(a,b,c)) for a,b,c in product(range(6),range(4),range(6))]
    assert {labels for labels,(_,i) in small if labels[0]==labels[2] and i<=15} == {
        (0,0,0),(0,1,0),(1,0,1)}
    assert max(F(idx,dim) for _,(dim,idx) in small if dim<=24) == F(21,20)

    data = run_path(str(root/"ch39-published-gauge-inputs.py"))
    scale = data["q"]("26.81725683017028")
    shared = data["model"]["target"](11431,scale)
    assert arb("3.29412973") < shared < arb("3.29412974")
    err = arb("9.08e-8")
    assert 15*35/(16*arb.pi()**2*(shared+err)) > arb("1.00925")
    assert 15*41/(16*arb.pi()**2*(shared+err)) > arb("1.18226")
    assert 15*21/(16*arb.pi()**2*(shared-err)) < arb(".60556")
    # If the weighted inverse couplings are only approximate, low
    # indices are still forced to the exact integer ratios by rational spacing.
    spread = 2*err/(shared-err)
    assert max(arb(16)/15,arb(2))*spread < arb(1)/(15*103)
    assert 6*103/(16*arb.pi()**2*(shared+err)) > arb("1.18226")
    print("primitive trace dimensions 35, 41 with left doublets; doubled trace dimension 56")
    print("direct simple-parent adjoint loop factor >1.00925, and >1.18226 with left doublets")


if __name__ == "__main__":
    verify()
