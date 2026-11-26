"""变化中心方向的长度，在平方弦距截止 224 处闭合最后一个边轨道。"""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

plane = run_path(str(Path(__file__).with_name("ch39-three-clock-plane-completion.py")))


def evaluate(t):
    points,distance,g0,delta,entry = plane["prepare"](t)
    def line(i,j,z):
        a,b,c,d = entry(i,j)
        return a+b*z,c+d*z
    def outer_coordinate(i,j):
        a,b,_,_ = entry(i,j)
        assert 0 not in b
        return -a/b
    x,w = outer_coordinate(4,2),outer_coordinate(8,0)
    def regular(j):
        a,b = line(4,j,x)
        return a+b*x
    A,B = regular(12),regular(15)
    assert 2*A>B and 2*A+B>0
    v = ((2*A-B)/3).sqrt()
    u = (-v-(2*A+B).sqrt())/2
    a,b = line(4,16,x)
    c0,c1,c2,c3 = entry(6,7)
    A2,A1,A0 = b*b-3*u*u*c3, 2*a*b-3*u*u*(c1+c2), a*a-3*u*u*c0
    disc = A1*A1-4*A2*A0
    assert disc>0 and 0 not in A2
    # 同一非退化二次式的较小实根在本区间连续。
    y = -A1/(2*A2)-disc.sqrt()/(2*abs(A2))
    r = (a+b*y)/(3*u)
    a,b = line(4,14,x)
    s = (a+b*w)/(3*u)
    a,b = line(6,9,y)
    last = a+b*w-3*r*s
    trace = arb(0)
    for i,p in enumerate(points):
        if sorted(p)==sorted(points[4]):
            z,squared = x,u*u+v*v+(u+v)**2
        elif sorted(p)==sorted(points[6]):
            z,squared = y,6*r*r
        elif sorted(p)==sorted(points[8]):
            z,squared = w,6*s*s
        else:
            z,squared = arb(0),arb(0)
        a,b,c,d = entry(i,i)
        trace += delta[i,i]+a+(b+c)*z+d*z*z+squared
    return last,trace/10**10


def verify():
    ctx.prec = 256
    points,distance,g0,_,base,slope = plane["shell_data"]()
    outer = plane["outer"]
    group = [[points.index(tuple(p[k] for k in perm)) for p in points]
             for perm in permutations(range(3))]
    for g in group:
        lookup = [outer.index(g[i]) for i in outer]
        assert all(m[i,j]==m[lookup[i],lookup[j]] for i in range(6) for j in range(6) for m in (base,slope))
    new_edges = {(i,j) for i,j in combinations(range(19),2) if 224<distance[i][j]<=244}
    assert new_edges == {tuple(sorted((g[8],g[0]))) for g in group}
    left,right = arb(1)/10,arb(1)/5
    assert evaluate(left)[0]>0>evaluate(right)[0]
    for _ in range(95):
        mid = (left+right)/2
        val = evaluate(mid)[0]
        assert 0 not in val
        if val>0:
            left = mid
        else:
            right = mid
    assert evaluate(left)[0]>0>evaluate(right)[0]
    box = left.union(right)
    error,trace = evaluate(box)
    assert error.contains(0) and trace>0
    inv_norm = max(sum(abs(x) for x in row) for row in g0.inv().tolist())
    assert trace*inv_norm < (1-arb(".098"))**2
    print("Cutoff 224, dimension at most 26, center parameter",box.str(21),
          "trace upper",trace.upper().str(12),flush=True)


if __name__ == "__main__":
    verify()
