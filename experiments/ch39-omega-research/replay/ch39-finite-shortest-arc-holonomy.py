"""Intrinsic shortest-arc transport, entirely over the rational numbers."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from runpy import run_path

rank = run_path(str(Path(__file__).with_name("ch39-finite-shell-transport-witness.py")))["rank"]


def eye(n):
    return [[int(i==j) for j in range(n)] for i in range(n)]


def mul(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def link(a,b):
    radius = sum(x*x for x in a)
    assert radius == sum(x*x for x in b)
    overlap = sum(x*y for x,y in zip(a,b))
    k = [[y*x-u*v for x,v in zip(a,b)] for y,u in zip(b,a)]
    k2 = mul(k,k)
    n = len(a)
    U = [[int(i==j)+Q(k[i][j],radius)+Q(k2[i][j],radius*(radius+overlap))
          for j in range(n)] for i in range(n)]
    assert mul(U,list(zip(*U)))==eye(n)
    assert [r[0] for r in mul(U,[[x] for x in a])]==list(b)
    return U


def loop(a,b,c):
    return mul(link(c,a),mul(link(b,c),link(a,b)))


def main():
    a,b,c,d = (0,1,3,4),(0,1,4,3),(0,3,1,4),(1,0,3,4)
    # The first six original vertices lie in the same three-coordinate hyperplane.
    # Embed one extra point with the same integer squared norm.
    paths = [(a,b,c),(a,c,d)]
    for path in paths:
        assert sorted(26-sum(x*y for x,y in zip(path[i],path[(i+1)%3])) for i in range(3)) == [1,4,7]
    h,k = [loop(*p) for p in paths]
    assert (sum(h[i][i] for i in range(4))-2)/2 == Q(6869,6885)
    # Tangent vectors at a have coordinates (x,-3y-4z,y,z).
    frame = [[1,0,0],[0,-3,-4],[0,1,0],[0,0,1]]
    def tangent(W):
        B = mul(W,frame)
        A = [B[0],B[2],B[3]]
        assert mul(frame,A)==B
        return A
    A,B = tangent(h),tangent(k)
    assert mul(A,B)!=mul(B,A)
    words = [eye(3)]
    level = [eye(3)]
    for length in range(1,5):
        level = [mul(W,V) for W,V in product(level,(A,B))]
        words.extend(level)
        dimension = rank([[v for row in W for v in row] for W in words])
        print("word length <=",length,"rational matrix span",dimension)
        if dimension==9:
            break
    assert dimension==9
    print("planar loop angle cosine exactly 6869/6885; two lifted loops do not commute")


if __name__=="__main__":
    main()
