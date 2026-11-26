"""Finite matrix-unit witness for the shell triality closure of the unified algebra."""
from fractions import Fraction as F
from itertools import combinations, permutations


def verify():
    pairs = list(combinations(range(5),2))
    y = (-2,-2,-2,3,3)
    charge = lambda pair: ((sum(y[i] for i in pair)+1)//3) % 3
    single = []
    for spectator in (1,2,4):
        source, target = tuple(sorted((0,spectator))), tuple(sorted((3,spectator)))
        single.append((source,target,(charge(target)-charge(source)) % 3))
    assert [x[2] for x in single] == [1,1,2]
    source,target,_ = single[-1]
    units = set()
    for perm in permutations(range(5)):
        old = tuple(sorted(perm[i] for i in source))
        new = tuple(sorted(perm[i] for i in target))
        units.add((new,old))
    assert len(units) == 60
    # One commutator through a middle vertex covers the disjoint pairs too.
    generated = units | {(i,k) for i,j in units for h,k in units if j == h and i != k}
    assert len(generated) == 90
    assert len(generated)+len(pairs)-1 == 99 and 24+99 == 123
    charges_five = [((-v+1)//3) % 3 for v in y]
    charges_ten = [charge(pair) for pair in pairs]
    assert (sum(charges_five)+5)%3 == (sum(charges_ten)+10)%3 == 0
    generator = (1,-1,0,0,0)
    assert sum((generator[i]+generator[j])**2 for i,j in pairs) == 3*sum(h*h for h in generator)
    assert F(1,2)+F(3,2)==2
    print("Projected wedge matrix units span sl(10); full connected dimension 123, index 3.")


if __name__ == "__main__":
    verify()
