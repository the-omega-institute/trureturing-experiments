"""Shortest-edge phase cost fixed by two Yukawa matrices on the relative triplet."""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from sympy import I, Rational as Q, symbols

folder = Path(__file__).resolve().parent
model = run_path(str(folder/"ch39-finite-relative-yukawa.py"))
sq = model["norm_squared"]


def main():
    # A shortest-edge square is K_{2,2}, with sides joining opposite parities.
    example = [(1+I*t)/(1-I*t) for t in map(Q,(0,1,2,3))]
    plus,minus = example[:2],example[2:]
    energy = sum(sq(p-m) for p in plus for m in minus)
    T,S = sum(plus)+sum(minus),sum(plus)-sum(minus)
    assert energy == 8-(sq(T)-sq(S))/2 and energy >= sq(S)
    pairs = list(combinations(range(4),2))
    points = list(permutations((0,1,3,4)))
    keys = lambda p:(tuple(sorted((p.index(0),p.index(1)))),(-1)**sum(p[i]>p[j] for i,j in pairs))
    short = {frozenset((p,q)) for p,q in combinations(points,2)
             if sum((x-y)**2 for x,y in zip(p,q)) == 2}
    assert short == {frozenset((p,q)) for p,q in combinations(points,2)
                     if keys(p)[0]==keys(q)[0] and keys(p)[1]!=keys(q)[1]}
    _,up,down = model["benchmark"]()
    def split(Y):
        even,odd = [],[]
        for p in pairs:
            c = tuple(sorted(set(range(4))-set(p)))
            even.append(-(Y[p]+Y[c])/2)
            odd.append(-(Y[p]-Y[c])/2)
        return sum(map(sq,even)),sum(map(sq,odd))
    u0,u1 = split(up)
    d0,d1 = split(down)
    assert u1 == 0 and d1>0
    scale = Q(3072,25)
    E0,penalty = 48-scale*(u0-d0),scale*d1
    E1,Ecent = E0+penalty,E0+10000*penalty
    bracket = lambda n:(int(n*10**8),int(n*10**8)+1)
    floor_all = 48-Q(3072,25)*Q(3,8)*sum(map(Q,(".00000615",".00311",".8616")))**2
    assert floor_all>Q("13.5444") and floor_all/24>Q(3,4)**2
    print("E at beta 1 and .01, lower/upper numerators over 1e8",bracket(E1),bracket(Ecent))
    print("Penalty coefficient",bracket(penalty),"/1e8; every single-cover realization at g=5 has a shortest chord above 3/4.")


if __name__=="__main__":
    main()
