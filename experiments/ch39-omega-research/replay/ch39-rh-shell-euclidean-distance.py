"""Optimal Euclidean two-point fit to the nine shell readings, already attained in a plane."""
from fractions import Fraction as Q
from itertools import combinations,product
from pathlib import Path
from runpy import run_path

distance=run_path(str(Path(__file__).with_name("ch39-rh-shell-metric-distance.py")))
roots=distance["roots"]
numbers=run_path(str(Path(__file__).with_name("ch39-rh-clock-integer-replay.py")))
I=numbers["Bounds"]


def main():
    lo,hi=Q("0.02817206520771579872"),Q("0.02817206520771579873")
    def error(e):
        short=distance["limits"](33,e)[1]
        medium=distance["limits"](25,e)[1]
        long=distance["limits"](18,e)[0]
        return I(*long)-I(*short)-I(*medium)
    assert error(lo).lo>0 and error(hi).hi<0
    def coordinate(n):
        return I(roots((1/(Q(n,34)-lo)-1)/2)[0],roots((1/(Q(n,34)-hi)-1)/2)[1])
    h,k=coordinate(33),coordinate(25)
    center=(I(2)*h+k)/I(3)
    mapping={0:I(0)-center,1:h-center,4:h+k-center}
    pts=[p for p in product(range(5),repeat=3) if sum(a*a for a in p)==17]
    embeddings=[]
    for p in pts:
        if len(set(p))==3:
            embeddings.append([mapping[i] for i in p])
        else:
            embeddings.append([I(Q(7,50)*(int(v==3)-Q(1,3))) for v in p])
    def sq(v):
        ends=[v.lo**2,v.hi**2]
        return I(0 if v.lo<=0<=v.hi else min(ends),max(ends))
    for (p,x),(q,y) in combinations(zip(pts,embeddings),2):
        n=17+sum(a*b for a,b in zip(p,q))
        if len(set(p))==len(set(q))==3 and n in (18,25,33):
            continue  # The defining transposition equality gives these three distances.
        d2=sum((sq(a-b) for a,b in zip(x,y)),I(0))
        measured=I(1)/(I(1)+d2)
        assert Q(n,34)-lo<measured.lo<=measured.hi<Q(n,34)+lo
    print("Exact Euclidean minimax error lies in",str(float(lo)),"certified between",lo,hi)
    print("Two-dimensional single-pole witness. Three transposition distances fix its lower bound.")


if __name__=="__main__":
    main()
