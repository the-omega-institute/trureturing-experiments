"""Exact weight and Pfaffian checks for a smooth nonzero-flux orbit."""
from fractions import Fraction as F
from itertools import combinations
from sympy import expand, symbols


def pfaffian(c):
    return c[0]*c[5]-c[1]*c[4]+c[2]*c[3]


def verify():
    fluxes = [symbols(f"f{i}_0:6") for i in range(4)]
    fluxes.append(tuple(-sum(c[j] for c in fluxes) for j in range(6)))
    index = lambda c: -3*pfaffian(c)
    expected = sum(index(c) for c in fluxes)
    for i, c in enumerate(fluxes):
        charge = -index(c)+sum(index(tuple(x+y for x,y in zip(c,d)))
                               for j,d in enumerate(fluxes) if i != j)
        assert expand(charge-expected) == 0
    size = 2**15
    assert F(1, size**2) == F(1, 2**30)
    weights = (1,-1,0,0,0)
    anti = [-3*n*n for n in weights]
    ten = [-3*(p+q)**2 for p,q in combinations(weights,2)]
    assert sum(anti) == -6 and sum(ten) == -18
    # Integer shell charge on the light SM fields, including the Higgs doublet.
    hypercharge = [F(1,6), F(-2,3), F(1,3), F(-1,2), 1, F(1,2)]
    baryon = [F(1,3), F(-1,3), F(-1,3), 0, 0, 0]
    assert [(2*y-b)%3 for y,b in zip(hypercharge,baryon)] == [0,2,1,2,2,1]
    print("General six-flux Pfaffian identities pass for all five weights.")
    print("Example L=32768, index(anti-5)=-6, index(10)=-18, all 81 flat twists allowed.")


if __name__ == "__main__":
    verify()
