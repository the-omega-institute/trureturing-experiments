"""Square classes of finite spin transports on the cutoff-four integer permutation shell."""
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import isqrt
from pathlib import Path
from runpy import run_path

cut = run_path(str(Path(__file__).with_name("ch39-finite-permutation-arc-cut.py")))
arc = cut["arc"]


def conj(a):
    return (a[0],-a[1],-a[2],-a[3])


def times(a,b):
    s,x = a[0],a[1:]
    t,y = b[0],b[1:]
    return (s*t-sum(i*j for i,j in zip(x,y)),
            *(s*y[k]+t*x[k]+x[(k+1)%3]*y[(k+2)%3]-x[(k+2)%3]*y[(k+1)%3] for k in range(3)))


def norm2(a):
    return sum(x*x for x in a)


def rotor(a,b):
    unit = times(conj(b),a)
    return (26+unit[0],*unit[1:])


def main():
    start = cut["start"]
    edges = []
    for a,b in combinations(permutations(start),2):
        d = 26-sum(x*y for x,y in zip(a,b))
        if d>4:
            continue
        if d==2:
            assert sum(x!=y for x,y in zip(a,b))==4
        else:
            assert sum(x!=y for x,y in zip(a,b))==2
        z = rotor(a,b)
        squareclass = {1:663,2:26,4:39}[d]
        assert norm2(z)==52*(52-d) and isqrt(norm2(z)//squareclass)**2*squareclass==norm2(z)
        U = arc["link"](a,b)
        for e in ((0,1,0,0),(0,0,1,0),(0,0,0,1)):
            transported = times(b,times(times(z,e),conj(z)))
            lhs = [Q(v,norm2(z)) for v in transported]
            assert lhs == [row[0] for row in arc["mul"](U,[[x] for x in times(a,e)])]
        edges.append(d)
    assert len(edges)==48 and 663*39==17*39**2
    paths = [[(0,1),(1,3)]*3,[(0,1),(3,4),((0,1),(3,4))]]
    for path,squareclass in zip(paths,(17,26)):
        a = start
        W = (1,0,0,0)
        for action in path:
            b = a
            for i,j in (action if isinstance(action[0],tuple) else [action]):
                b = tuple(j if x==i else i if x==j else x for x in b)
            W = times(rotor(a,b),W)
            a = b
        assert a==start and W[0]!=0
        n = norm2(W)
        assert n == squareclass*isqrt(n//squareclass)**2
        print("closed path length",len(path),"spin trace squared",Q(4*W[0]**2,n),"squareclass",squareclass)


if __name__=="__main__":
    main()
