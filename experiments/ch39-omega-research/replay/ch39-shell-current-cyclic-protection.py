"""Local triangle words that protect an exchange against every cyclic-equivariant kernel."""
from collections import Counter
from itertools import combinations
from math import gcd
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
shell = run_path(str(root/'ch39-shell-current-word-floor.py'))
points, links = shell['points'], shell['transfer']
pole, radius2, dot = shell['pole'],shell['radius2'],shell['dot']


def primitive(pair):
    a,b = pair
    g = gcd(a,b)
    return (a//g,b//g)


def multiply(u,v):
    a,b = u
    c,d = v
    return primitive((a*c-3*b*d,a*d+b*c))


def power(u,n):
    result = (1,0)
    while n:
        if n%2: result = multiply(result,u)
        n //= 2
        if n: u = multiply(u,u)
    return result


def link(i,j):
    x,y = points[i],points[j]
    det = sum(pole[k]*(x[(k+1)%3]*y[(k+2)%3]-x[(k+2)%3]*y[(k+1)%3]) for k in range(3))
    raw = (radius2*(radius2+dot(pole,x)+dot(x,y)+dot(y,pole)),-17*det)
    if 2*(radius2-dot(x,y))==200: raw = tuple(-t for t in raw)
    return power(primitive(raw),5)


def classify():
    words = {edge:link(*edge) for edge in links}
    orbit = [min(v,v[1:]+v[:1],v[2:]+v[:2]) for v in points]
    counts = Counter(orbit)
    assert sorted(counts.values())==[1]+[3]*6
    identity = nonlocal_only = 0
    local_budget = Counter()
    transport_budget = Counter()
    for i,j,k in combinations(range(len(points)),3):
        if not all(e in words for e in ((i,j),(j,k),(k,i))): continue
        c = multiply(words[i,j],multiply(words[j,k],words[k,i]))
        if c==(1,0):
            identity += 1
            continue
        assert power(c,12)!=(1,0)
        # c squared is in Q(sqrt(-3)). Its being torsion would force c^12=1.
        directed = (c[0],abs(c[1]))
        assert directed[1]>0
        winding = directed
        turns = 1
        while winding[1]>0:
            winding = multiply(winding,directed)
            turns += 1
        assert winding[1]<0
        transport_budget[turns]+=1
        if len({orbit[i],orbit[j],orbit[k]})==3:
            nonlocal_only += 1
        else:
            local_budget[turns]+=1
    assert local_budget=={2:61,6:3,7:24,8:3,9:6,13:6,15:18,18:9,23:12,58:6,68:6}
    assert sum(local_budget.values())==154 and nonlocal_only==372 and identity==3
    assert sum(transport_budget.values())==526 and max(transport_budget)==916
    print('All equivariant kernels, local eligible 154, local ineligible 375.')
    print('Local minimum turn counts and multiplicities:',sorted(local_budget.items()))
    print('Transport: 526 eligible up to 916 turns; three identity triangles.')


if __name__=='__main__': classify()
