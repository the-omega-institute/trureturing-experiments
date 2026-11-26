"""Center charge on paired diagonal embeddings and finite prime-root avoidance."""
from fractions import Fraction as F
from itertools import product
from math import ceil

import sympy as sp


def center(f,g,a,b):
    n = len(f)
    charges = (a,-a,b,-b,0)
    conjugate = sum(f[-x%n] for x in charges)
    exterior = sum(g[(x+y)%n] for i,x in enumerate(charges) for y in charges[i+1:])
    return -conjugate+2*exterior


def even_tables(p):
    for half in product(range(5),repeat=p//2):
        vals = (0,*half)
        yield tuple(vals[min(i,p-i)] for i in range(p))


def identities():
    for p in (3,7):
        for g in even_tables(p):
            rows = [center(g,g,a,1) for a in range(p)]
            assert (sum(rows)-2*p*g[1])%5==0
            assert not g[1] or any(z%5 for z in rows)
            for f in (g,tuple(0 for _ in g)):
                assert all((center(f,g,a,0)-2*(g[a]-f[a]))%5==0 for a in range(p))
                for a in range(p):
                    b=1
                    half = -f[a]-f[b]+2*(g[(a+b)%p]+g[(a-b)%p]+g[a]+g[b])
                    assert center(f,g,a,b)==2*half
    z=sp.symbols("z")
    for p,q in ((7,11),(11,13),(13,17)):
        assert sp.gcd(sp.cyclotomic_poly(p,z),sp.cyclotomic_poly(q,z))==1


def tent(delta,C):
    r = ceil(16*C/delta)
    e = delta/(24*C*r)
    assert C*(F(2,r+1)+3*r*e)<delta/4


if __name__ == "__main__":
    identities()
    tent(F(1,50),F(1234567))
    print("Center identities, disjoint primitive prime roots, and the finite homotopy-gap constants verified.")
