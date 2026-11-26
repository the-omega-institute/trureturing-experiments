"""Permutation-covariant flavor links on every shortest core edge, including conjugate charge."""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx


folder=Path(__file__).resolve().parent
calc=run_path(str(folder/"ch39-shell-gw-rh-flavor.py"))
ctx.prec=240
model=calc["model"]
points,faces=model["data"]["maps"]()


def link(v,w,p,q):
    unit=[[acb(int(i==j)) for j in range(3)] for i in range(3)]
    if 26-sum(x*y for x,y in zip(v,w))!=1:
        return unit
    i,j=v.index(3),v.index(4)
    unit[i][i]=acb(p.cos(),p.sin())*q.cos()
    unit[j][j]=unit[i][i].conjugate()
    unit[i][j]=-q.sin()
    unit[j][i]=q.sin()
    return unit


def equivariant_matrix(base,p,q):
    # The block in row x, column y transports from y to x.
    t={(i,j):link(points[j],points[i],p,q) for i in range(6) for j in range(6)}
    return [[base[i//3][j//3]*t[(i//3)%6,(j//3)%6][i%3][j%3]
             + (arb(3)/4 if i==j and i<18 else -arb(3)/4 if i==j else 0)
             for j in range(36)] for i in range(36)]


if __name__=="__main__":
    p,q=arb(3)/5,arb(4)/7
    for v,w in combinations(points.values(),2):
        uv=link(v,w,p,q)
        vu=link(w,v,p,q)
        assert all(acb(uv[i][j]).overlaps(acb(vu[j][i]).conjugate()) for i in range(3) for j in range(3))
        for perm in permutations(range(3)):
            observed=link(tuple(v[i] for i in perm),tuple(w[i] for i in perm),p,q)
            assert all(acb(observed[i][j]).overlaps(acb(uv[perm[i]][perm[j]])) for i in range(3) for j in range(3))
    angles=[(1-arb(d)/26).acos() for d in (3,5,7)]
    f,g=angles[1]-angles[0],angles[2]-angles[0]
    examples=(1,4,13,25,255,719)
    for n in examples:
        h=acb.zeta_zero(n).imag
        result=[]
        for sign in (1,-1,0):
            base=model["make"](sign,True)
            matrix=equivariant_matrix(base,(-1 if sign<0 else 1)*h*f,h*g)
            index=18-model["count"](matrix,arb(0))
            result.append(index)
        print(n,result)
