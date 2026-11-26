"""Shortest affine heap words for antipodes of the nine positive Behrend points."""
from collections import Counter
from itertools import combinations, permutations, product
from fractions import Fraction as Q

outer=(0,1,4)
inner=(2,2,3)
positive=[outer]*3+[(0,4,1)]+[(4,1,0)]*3
negative=[inner]*3+[(3,2,2)]*2


def heap(a,b,c):
    return tuple(x-y+z for x,y,z in zip(a,b,c))


def reflect(p):
    # Translate the single doubled-center relation so the point belongs to it.
    order=next(order for order in permutations(range(3)) if
               tuple((outer if sum(p)==5 else inner)[j] for j in order)==p)
    coeff=Counter(tuple(x[j] for j in order) for x in positive)
    coeff.subtract(tuple(x[j] for j in order) for x in negative)
    coeff[p]-=1
    plus=list(coeff.elements())
    minus=list((-coeff).elements())
    assert len(plus)==len(minus)+1
    t=plus.pop()
    for x,y in zip(minus,plus):
        t=heap(t,x,y)
    return t,sum(abs(x) for x in coeff.values())


def main():
    shell=[v for v in product(range(5),repeat=3) if sum(x*x for x in v)==17]
    assert all(sum(x[i] for x in positive)==sum(x[i] for x in negative) for i in range(3))
    for p in shell:
        q,n=reflect(p)
        assert q==tuple(-x for x in p) and n==(11 if sum(p)==5 else 13)
    for p,q in combinations(shell,2):
        reflected,_=reflect(q)
        add=sum((x-y)**2 for x,y in zip(p,reflected))
        sub=sum((x-y)**2 for x,y in zip(p,q))
        assert add+sub==68 and Q(add,add+sub)==Q(17+sum(x*y for x,y in zip(p,q)),34)
    print("All 36 two-channel measurements are affine and integral.")
    print("Antipodes need 11 or 13 leaf occurrences, respectively; the center has no affine integer word.")


if __name__=="__main__":
    main()
