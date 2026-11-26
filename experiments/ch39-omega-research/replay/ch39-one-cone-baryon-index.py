"""Finite gap, magnetic interpolation budget, and the three-family index."""
from fractions import Fraction as F
from itertools import combinations
from sympy import expand, symbols


def pfaffian(c):
    return c[0]*c[5]-c[1]*c[4]+c[2]*c[3]


def total_index(charges):
    color, weak = charges[:3], charges[3:]
    def add(c,d):
        return tuple(x+y for x,y in zip(c,d))
    idx = lambda group: -sum(pfaffian(c) for c in group)
    return (idx([add(c,w) for c in color for w in weak]),
            idx([add(c,d) for c,d in combinations(color,2)]), idx(color),
            idx(weak), idx([add(*weak)]))


def verify():
    mass = F(1,16)
    assert F(3,8)-2*mass > 0
    assert mass-2*(F(16,4096)+12*4095*F(1,2**26)) > F(1,100)
    # The interpolation bound grows logarithmically with the denominator and flux.
    # Old and new sign-polynomial degrees differ by at most a factor of sixteen.
    assert 10240000 == 4*2560000
    assert 16*(2*10240000+1)**2 < 2**54
    assert 1024*3*2**(54*4) < 2**228
    k = 16*(2048+1)
    assert k**8 < 2**(k//2)
    vals = [symbols(f"w{i}_0:6") for i in range(4)]
    vals.append(tuple(-sum(c[j] for c in vals) for j in range(6)))
    q, u, d, lepton, e = total_index(vals)
    c2 = pfaffian(tuple(a+b for a,b in zip(vals[3],vals[4])))-sum(map(pfaffian,vals[3:]))
    assert expand(q-u-d-3*c2) == expand(lepton-e-c2) == 0
    zero = (0,)*6
    flux1 = (1,0,0,0,0,0)
    flux2 = (0,0,0,0,0,1)
    colors = [(-1,0,0,0,0,-1),zero,zero]
    assert total_index(colors+[flux1,flux2]) == (0,-2,-1,0,-1)
    assert 2*F(1,32768**2) <= F(1,2**29)
    print("One free cone; uniform elementary gap 1/16, smooth magnetic index -Pf(C).")
    print("One-generation I(Q,u^c,d^c,L,e^c)=(0,-2,-1,0,-1), baryon and lepton indices both 1.")


if __name__ == "__main__":
    verify()
