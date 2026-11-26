"""Finite witness of the traceless 3+2 block gauge algebra on two shell circuits."""
from fractions import Fraction as Q
from math import gcd, isqrt
from pathlib import Path
from runpy import run_path

quat = run_path(str(Path(__file__).with_name("ch39-finite-quaternion-loop-field.py")))
cut = quat["cut"]
p = 101
imag = 10
sqrt17 = next(i for i in range(p) if i*i%p==17)


def reduce(x):
    v = Q(x)
    return v.numerator*pow(v.denominator,-1,p)%p


def pair(path):
    a = cut["start"]
    z = (1,0,0,0)
    for i,j in path:
        b = tuple(j if x==i else i if x==j else x for x in a)
        z = quat["times"](quat["rotor"](a,b),z)
        a = b
    assert a == cut["start"]
    z = tuple(t//gcd(*z) for t in z)
    n = quat["norm2"](z)
    if isqrt(n)**2==n:
        scale = isqrt(n)
    else:
        assert 17*isqrt(n//17)**2==n
        scale = sqrt17*isqrt(n//17)
    w,x,y,t = [(u*pow(scale,-1,p))%p for u in z]
    spin = [[w+imag*x,y+imag*t],[-y+imag*t,w-imag*x]]
    vector = cut["walk"](path)
    return [[reduce(vector[i][j]) if max(i,j)<3 else spin[i-3][j-3]%p if min(i,j)>=3 else 0
             for j in range(5)] for i in range(5)]


def expand(gens):
    basis = {}
    def add(row):
        v = [x%p for x in row]
        for i in sorted(basis):
            v = [(a-v[i]*b)%p for a,b in zip(v,basis[i])]
        if any(v):
            i = next(j for j in range(len(v)) if v[j])
            inv = pow(v[i],-1,p)
            basis[i] = [x*inv%p for x in v]
            return True
        return False
    flat = lambda v:[x for r in v for x in r]
    for g in gens:
        assert add(flat(g))
    active = list(gens)
    sizes = [len(active)]
    for _ in range(10):
        fresh = []
        for K in active:
            for G in gens:
                L,R = cut["mul"](K,G),cut["mul"](G,K)
                Z = [[(a-b)%p for a,b in zip(l,r)] for l,r in zip(L,R)]
                if add(flat(Z)):
                    fresh.append(Z)
        if not fresh:
            return sizes
        sizes.append(len(basis))
        active = fresh
    raise AssertionError("unexpected bracket depth")


def main():
    square = [(0,1),(3,4)]*2
    hexagon = [(0,1),(1,3)]*3
    A,Ai,B,Bi = [pair(s) for s in (square,square[::-1],hexagon,hexagon[::-1])]
    for v,z in [(A,Ai),(B,Bi)]:
        assert [[x%p for x in row] for row in cut["mul"](v,z)] == cut["eye"](5)
    generators = [[[x-y for x,y in zip(rl,rr)] for rl,rr in zip(v,z)] for v,z in [(A,Ai),(B,Bi)]]
    assert expand(generators)[-1] == 6
    even = [[imag*(x+y) for x,y in zip(rl,rr)] for rl,rr in zip(A,Ai)]
    total = sum(even[i][i] for i in range(5))*pow(5,-1,p)
    for i in range(5):
        even[i][i] -= total
    dims = expand(generators+[even])
    assert dims[-1]==12
    print("Paired tangent and spin block bracket ranks modulo 101",dims)
    # Primitive central charges 2 and -3. Usual choice of sign and sixth-unit.
    Y = [-Q(1,3)]*3+[Q(1,2)]*2
    assert sum(Y)==0 and sum(t*t for t in Y)==Q(5,6)
    def trace_ratio(triplet, doublet):
        weak = doublet/2
        hypercharge = triplet/3 + doublet/2
        return weak/(weak+hypercharge)
    assert trace_ratio(Q(1),Q(1))==Q(3,8)
    assert trace_ratio(Q(2),Q(1))==Q(3,10)
    assert trace_ratio(Q(1),Q(2))==Q(3,7)


if __name__=="__main__":
    main()
