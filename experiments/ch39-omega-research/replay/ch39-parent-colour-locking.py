"""Minimum scalar index for locking the two optimal SU(21) colour embeddings."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations, product
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent


def splitting():
    block = run_path(str(root/"ch39-simple-parent-indices.py"))["block"]
    # Each defining space has at most 21 dimensions; conjugation is immaterial.
    choices = []
    for a,b,c,l,r in product(range(4),range(3),range(4),range(4),range(4)):
        labels = (a,b,c,l,r)
        d,index = block(labels)
        if c<=a and d<=21 and any(index):
            choices.append((labels,d,index))
    assert block((4,0,0,0,0))[0]>21 and block((0,3,0,0,0))[0]>21
    assert block((0,0,0,4,0))[1][1]>16

    @lru_cache(None)
    def fill(j,space,targets):
        if not any(targets):
            return ((),)
        results = []
        for k in range(j,len(choices)):
            w,d,index = choices[k]
            rem = tuple(x-y for x,y in zip(targets,index))
            if space>=d and min(rem)>=0:
                results.extend((w,)+s for s in fill(k,space-d,rem))
        return tuple(results)

    groups = set()
    for parts in fill(0,42,(15,16,14)):
        if not any(w[3]%2 for w in parts):
            continue
        for side in product(range(2),repeat=len(parts)):
            split = [tuple(w for w,j in zip(parts,side) if j==i) for i in range(2)]
            if all(sum(block(w)[0] for w in group)<=21 for group in split):
                groups.add(tuple(sorted(split)))
    witness = (((0,0,0,0,1),)*2+((0,0,0,2,2),(1,0,0,1,0)),((1,1,0,0,0),))
    assert groups == {witness}


def small_parent_representations():
    # Weyl dimension and twice the Dynkin index for SU(21); columns are
    # heights of Young columns, with no determinant column.
    n = 21
    denominator = prod(j-i for i,j in combinations(range(n),2))

    def representation(columns):
        rows = [sum(h>i for h in columns) for i in range(n)]
        dim = prod(rows[i]-rows[j]+j-i for i,j in combinations(range(n),2))//denominator
        boxes = sum(rows)
        numerator = (n*(n*boxes+sum(z*(z+1-2*i) for i,z in enumerate(rows,1)))-boxes**2)*dim
        divisor = n*(n*n-1)
        assert numerator%divisor == 0
        return dim,F(numerator//divisor,2)

    def grow(columns):
        for h in range(columns[-1] if columns else 1,n):
            w = (*columns,h)
            dim,idx = representation(w)
            if idx<=115:
                yield w,dim,idx
                yield from grow(w)
    small = {w:(d,i) for w,d,i in grow(())}
    canonical = lambda w:min(w,tuple(n-h for h in reversed(w)))
    entries = {canonical(w):v for w,v in small.items()}
    assert entries == {(1,):(21,F(1,2)), (2,):(210,F(19,2)), (1,1):(231,F(23,2)),
                       (1,20):(440,F(21)), (3,):(1330,F(171,2))}
    # A nontrivial self-conjugate representation has dimension and index
    # at least the adjoint's, since its columns contain a complementary pair.
    for h in range(1,11):
        dim,idx = representation((h,n-h))
        assert dim>=440 and idx>=21
    assert 440*21/2 > 2415
    return entries


def colour_squares():
    # Character of the (1,1,0) twenty = 4 x wedge^2(4) - wedge^3(4).
    unit = [tuple(int(i==j) for j in range(4)) for i in range(4)]
    def add(*ws):
        return tuple(map(sum,zip(*ws)))
    weights = Counter(add(a,b,c) for a in unit for b,c in combinations(unit,2))
    weights.subtract(add(*bs) for bs in combinations(unit,3))
    assert min(weights.values())>=0
    modes = list(weights.elements())
    assert len(modes)==20
    rho = (3,2,1,0)
    counts = []
    for symmetric in (False,True):
        pairs = Counter(add(a,b) for a,b in combinations(modes,2))
        if symmetric:
            pairs.update(add(a,a) for a in modes)
        multiplicities = []
        for highest in ((3,1,1,1),(2,2,2,0),(2,2,1,1)):
            value = 0
            for t in permutations(rho):
                sign = (-1)**sum(t[i]<t[j] for i,j in combinations(range(4),2))
                value += sign*pairs[tuple(a+b-c for a,b,c in zip(highest,rho,t))]
            multiplicities.append(value)
        counts.append(multiplicities)
    assert counts == [[0,0,1],[1,1,0]]  # ten, conjugate ten, six


def verify():
    splitting()
    entries = small_parent_representations()
    costs = {}
    for (w,(d,i)),(v,(e,j)) in product(entries.items(),repeat=2):
        cost = max(e*i,d*j)
        if cost<=2415:
            costs[tuple(sorted((w,v)))] = cost
    assert len(costs)==7
    assert set(costs.values()) == set(map(F,("21/2","399/2","483/2","441","3591/2","1995","2415")))
    colour_squares()
    assert F(-11*21+2415,3)==728
    # At this threshold use the target and allowance from the single-word source.
    data = run_path(str(root/"ch39-published-gauge-inputs.py"))
    inverse = (data["model"]["target"](11431,data["q"]("26.81725683017028"))+arb("9.08e-8"))/15
    assert (8*arb.pi()**2*inverse/728).exp() < arb("1.02411")
    print("unique optimal split; minimum complex scalar index 2415; one-loop coefficient >=728")


if __name__ == "__main__":
    verify()
