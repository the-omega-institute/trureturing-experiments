#!/usr/bin/env python3
"""Thankachy, Chandran, Tuite, Thomas, Di Stefano, Erskine, "On the Vertex Position Number of Graphs",
arXiv:2209.00359v2; Discuss. Math. Graph Theory 44 (2024) 1169-1188. Problem 2.4: "Is the ratio vp(G)/vp^-(G)
bounded for connected graphs?"

Family G_m (m >= 2): N = 4m, a_j = 2m+2j, b_j = 2m+j+1 (1 <= j <= m). Paths from a common vertex y:
H_0..H_N, P_{i,0..N} (1<=i<=m), Q_{j,0..a_j} (1<=j<=m), all index-0 vertices equal to y. Extra edges
H_{b_j} - Q_{j,b_j} and Q_{j,a_j} - P_{i,a_j} for all i, j.

Checks for m = 2..14 (exact BFS):
 (1) vertex and edge counts 7m^2+5m+1 and 8m^2+6m, connectivity;
 (2) every vertex lies on one of the 2m+1 paths, and each path is geodesic from y (dist = index), so
     p_y <= 2m+1 (an y-position set meets each geodesic path at most once, and y forces a singleton);
     the 2m+1 neighbours of y form a y-position set, so p_y = 2m+1;
 (3) the set S = {P_{i,a_j+1} : 1<=i<=m, 1<=j<=m-1} lies on one distance sphere around x = H_N (distance 2m+1),
     hence is an x-position set of size m(m-1);
 (4) by brute force for m = 2, 3 the position condition of S and of the neighbours of y.
Ratio: vp/vp^- >= m(m-1)/(2m+1), unbounded; m = 14 gives > 6."""
import sys, itertools
from collections import deque
bad=0
def rep(l,ok):
    global bad; print(l,ok); bad+=0 if ok else 1
def build(m):
    N=4*m; a={j:2*m+2*j for j in range(1,m+1)}; b={j:2*m+j+1 for j in range(1,m+1)}
    V={}; paths=[]
    def v(name,k):
        if k==0: name='y'
        key=(name,k if name!='y' else 0)
        if key not in V: V[key]=len(V)
        return V[key]
    E=set()
    def add(u,w): E.add((min(u,w),max(u,w)))
    specs=[('H',N)]+[(('P',i),N) for i in range(1,m+1)]+[(('Q',j),a[j]) for j in range(1,m+1)]
    for name,L in specs:
        p=[v(name,k) for k in range(L+1)]; paths.append(p)
        for k in range(L): add(p[k],p[k+1])
    for j in range(1,m+1):
        add(v('H',b[j]),v(('Q',j),b[j]))
        for i in range(1,m+1): add(v(('Q',j),a[j]),v(('P',i),a[j]))
    adj=[[] for _ in V]
    for u,w in E: adj[u].append(w); adj[w].append(u)
    return V,E,adj,paths,v,a
def bfs(adj,s):
    d=[-1]*len(adj); d[s]=0; q=deque([s])
    while q:
        u=q.popleft()
        for w in adj[u]:
            if d[w]<0: d[w]=d[u]+1; q.append(w)
    return d
def is_position(adj,x,S):
    dx=bfs(adj,x); S=list(S)
    for u in S:
        du=bfs(adj,u)
        for w in S:
            if u!=w and dx[u]+du[w]==dx[w]: return False
    return True
for m in range(2,15):
    V,E,adj,paths,v,a=build(m)
    ok=len(V)==7*m*m+5*m+1 and len(E)==8*m*m+6*m
    y=v('y',0); dy=bfs(adj,y); ok&=min(dy)>=0
    cover=set().union(*map(set,paths))==set(range(len(V)))
    geo=all(dy[p[k]]==k for p in paths for k in range(len(p)))
    x=v('H',4*m); dx=bfs(adj,x)
    S=[v(('P',i),a[j]+1) for i in range(1,m+1) for j in range(1,m)]
    sphere=len(set(S))==m*(m-1) and all(dx[s]==2*m+1 for s in S)
    ok&=cover and geo and sphere
    if m<=3:
        ok&=is_position(adj,x,S) and is_position(adj,y,adj[y]) and len(adj[y])==2*m+1
    rep(f'm={m}: n={len(V)} ratio>={m*(m-1)}/{2*m+1}={m*(m-1)/(2*m+1):.3f}',ok)
print('bad=',bad); sys.exit(1 if bad else 0)
