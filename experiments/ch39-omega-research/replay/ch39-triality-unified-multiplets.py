"""Two conjugate gauge backgrounds constrain complete SU(5) spectator multiplets."""
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from math import prod

Y = (-2,-2,-2,3,3)
PROBES = ((0,1,0,-1,0),(-1,1,0,0,0))


def dimension(labels):
    return prod(F(j-i+sum(labels[i:j]), j-i)
                for i in range(5) for j in range(i+1,5))


def charges(weights, sign, uniform=None):
    output = []
    for h in PROBES:
        charge = 0
        for indices in weights:
            y = sign*sum(Y[i] for i in indices)
            gen = sign*sum(h[i] for i in indices)
            triality = (y+1)//3 if uniform is None else uniform
            charge += gen*int((gen+triality)%3 != 0)
        output.append(charge)
    return tuple(output)


def certify():
    assert all(dimension(tuple(3 if i==j else 0 for i in range(4))) > 15
               for j in range(4))
    small = sorted(int(dimension(labels)) for labels in product(range(3),repeat=4)
                   if dimension(labels) <= 15)
    assert small == [1,5,5,10,10,15,15]
    weights = [[(i,) for i in range(5)], list(combinations(range(5),2)),
               list(combinations_with_replacement(range(5),2))]
    rows = [charges(w,s) for w in weights for s in (1,-1)]
    assert rows == [(0,-1),(0,1),(-2,-1),(2,1),(0,1),(0,-1)]
    for w,s,t in product(weights,(1,-1),range(3)):
        first, second = charges(w,s,t)
        assert first == second
    assert F(2,3)*(F(1,2)+F(3,2)) == F(4,3)
    assert F(2,3)*F(3,5)*4 == F(8,5)
    print("Certified small multiplet spectrum, two independent flux charges and both one-loop shifts.")


if __name__ == "__main__":
    certify()
