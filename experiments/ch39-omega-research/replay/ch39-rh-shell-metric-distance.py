"""Best nine-point metric fit for any positive normalized resolvent response."""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import isqrt


def roots(q):
    scale=2**100
    n=isqrt(q.numerator*scale**2//q.denominator)
    return Q(n,scale),Q(n+1,scale)


def limits(n,e):
    k=Q(n,34)
    return roots(1/(k+e)-1),roots(1/(k-e)-1)


def main():
    points=[p for p in product(range(5),repeat=3) if sum(i*i for i in p)==17]
    numer={(p,q):17+sum(i*j for i,j in zip(p,q)) for p,q in permutations(points,2)}
    arcs=set(numer.values())
    lo,hi=Q("0.023300655712046"),Q("0.023300655712049")
    bounds={n:limits(n,lo) for n in arcs}
    witness=((0,1,4),(2,2,3),(4,1,0))
    a,b,c=witness
    assert (numer[a,c],numer[a,b],numer[b,c])==(18,31,27)
    assert bounds[18][0][0]>bounds[31][1][1]+bounds[27][1][1]
    right={n:limits(n,hi) for n in arcs}
    assert right[18][0][1]<right[31][1][0]+right[27][1][0]
    assert 4*min(v[1][0] for v in bounds.values())>max(v[0][1] for v in bounds.values())
    ties=0
    for count in (2,3,4):
        for path in permutations(points,count):
            edges=[numer[p,q] for p,q in zip(path,path[1:])]
            n=numer[path[0],path[-1]]
            if n==18 and sorted(edges)==[27,31]:
                ties+=1
            else:
                assert sum(bounds[v][1][0] for v in edges)>bounds[n][0][1]
    assert ties==12
    print("Metric minimax readout error between",str(float(lo)),"and",str(float(hi)))
    print("Six unoriented tight paths; all other paths have positive slack.")


if __name__=="__main__":
    main()
