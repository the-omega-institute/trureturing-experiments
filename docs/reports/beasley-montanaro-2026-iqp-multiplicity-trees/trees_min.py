import sys
from collections import defaultdict
def canon(n, E):
    g=[[] for _ in range(n)]
    for u,v in E: g[u].append(v); g[v].append(u)
    def enc(x,p): return "("+"".join(sorted(enc(y,x) for y in g[x] if y!=p))+")"
    return min(enc(r,-1) for r in range(n))
def table(E):
    m=len(E); t=defaultdict(int)
    for s in range(1<<m):
        bd=0
        for e in range(m):
            if s>>e&1: u,v=E[e]; bd^=(1<<u)|(1<<v)
        t[(bin(s).count('1'),bin(bd).count('1'))]+=1
    return tuple(sorted(t.items()))
level={canon(2,[(0,1)]):[(0,1)]}
for n in range(3,14):
    nxt={}
    for E in level.values():
        for v in range(n-1):
            F=E+[(v,n-1)]; c=canon(n,F)
            if c not in nxt: nxt[c]=F
    level=nxt
    groups=defaultdict(list)
    for E in level.values(): groups[table(E)].append(E)
    coll=[g for g in groups.values() if len(g)>1]
    print(n,'trees',len(level),'colliding classes',len(coll),'trees in collisions',sum(len(g) for g in coll),flush=True)
    if coll and n<=11:
        for g in coll[:3]: print('   ',g)
    if n>=12: break
