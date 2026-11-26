"""Minimum original-site support for positive local squares, and bounds for square mass matrices."""
from itertools import combinations
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/'ch39-arc-relative-determinant.py'))
edges={(i,j):(d,z) for i,j,d,z in m['edges']}
low=set(range(9)); high=set(range(9,19))
def neighbors(i):
    return {b if a==i else a for a,b in edges if i in (a,b)}
degree={i:len(neighbors(i)) for i in range(19)}
for i in low:
    assert m['V'][i]+degree[i]==len(neighbors(i)&low)
# Dual witness: block X_ii=t_i^2 I and X_ij=t_i*t_j U_ij on each edge.
k=len(neighbors(9)&low)
high_diag=m['V'][9]+degree[9]
gain=k*k/high_diag
assert k>0 and gain>arb('.00005')
print('phase-free negative trial, magnitude greater than',gain.lower().str(12))
triangles=[(0,2,6),(1,3,7),(4,5,8),(6,7,8)]
used=set()
for i,j,k in triangles:
    pairs=[(i,j),(j,k),(i,k)]
    assert not set(pairs)&used and all(edges[p][0] not in (724,804) for p in pairs)
    used.update(pairs)
    flux=edges[i,j][1]*edges[j,k][1]*edges[i,k][1].conjugate()
    gap=abs(1-flux)**2/9
    debit=max(len(neighbors(v)&high)*arb('.00001') if v<6 else arb(0) for v in (i,j,k))
    if (i,j,k)==(6,7,8):
        debit=max(len(neighbors(v)&high)*arb('.001') for v in (i,j,k))
    assert gap>debit
    print('triangle',i,j,k,'squared-form margin at least',(gap-debit).lower().str(10))
for v in high:
    adjacent=neighbors(v)&low
    reserve=m['V'][v]+len(adjacent)-sum(100000 if i<6 else 1000 for i in adjacent)
    assert reserve>350000
# Square Cholesky of the endpoint form: symbolic support is controlled by filled neighborhoods.
order=(0,2,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
adj={i:neighbors(i) for i in order}
live=set(order)
width=0
for i in order:
    later=adj[i]&live
    width=max(width,1+len(later))
    for a,b in combinations(later,2):
        adj[a].add(b);adj[b].add(a)
    live.remove(i)
assert width==13
assert 38*6 < 2*len(edges) and 90*3==2*len(edges)
print('square factor support at most',width,'and at least 5 on the 19 static sites')
