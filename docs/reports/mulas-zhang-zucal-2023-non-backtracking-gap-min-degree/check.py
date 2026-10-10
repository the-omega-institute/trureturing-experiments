#!/usr/bin/env python3
"""Mulas, Zhang, Zucal, arXiv:2303.00373v2, Section 5, Conjecture 5.3: for a simple graph G with minimum degree
delta >= 2 that is not the cycle graph, the spectral gap from 1 of the non-backtracking Laplacian,
eps = min |lambda| over the spectrum of the non-backtracking transition matrix D^{-1}A, satisfies eps <= 1/(delta-1).

Checks (floating point, numpy): wheels W_N (N = 4..11; the source reports eps = 1/2 at N = 4, then decreasing),
complete graphs K_4..K_8, and 3000 random graphs on 4..12 vertices with minimum degree >= 2. For each graph it
checks eps <= 1/(delta-1) and, when delta >= 3, the stronger statement that the transition matrix has a positive
real eigenvalue in (0, 1/(delta-1)] (the proof's symmetric-pencil crossing)."""
import numpy as np, itertools, random, sys
random.seed(7); np.random.seed(7)
def P_of(n, edges):
    adj={v:set() for v in range(n)}
    for u,v in edges: adj[u].add(v); adj[v].add(u)
    arcs=[(u,v) for u in range(n) for v in sorted(adj[u])]
    idx={a:i for i,a in enumerate(arcs)}
    P=np.zeros((len(arcs),len(arcs)))
    for (u,v) in arcs:
        d=len(adj[v])
        for w in adj[v]:
            if w!=u: P[idx[(u,v)], idx[(v,w)]] = 1.0/(d-1)
    return P, adj
bad=0; tested=0
def test(n, edges, label=None):
    global bad, tested
    P,adj=P_of(n,edges); degs=[len(adj[v]) for v in range(n)]
    d=min(degs)
    if d<2: return
    ev=np.linalg.eigvals(P); eps=min(abs(ev))
    pos=[e.real for e in ev if abs(e.imag)<1e-9 and 0<e.real<=1/(d-1)+1e-9]
    ok = eps <= 1/(d-1)+1e-9 and (d==2 or len(pos)>0)
    tested+=1
    if not ok: bad+=1; print('FAIL',label,n,edges,eps,d)
    return eps,d,pos
# wheels (source example), complete graphs, random graphs with delta>=3 incl. irregular
for N in range(4,12):
    hub=N-1; e=[(i,(i+1)%(N-1)) for i in range(N-1)]+[(i,hub) for i in range(N-1)]
    r=test(N,e,'wheel'); print('wheel',N,'eps=%.4f'%r[0],'delta',r[1])
for n in range(4,9): test(n,list(itertools.combinations(range(n),2)),'K')
cnt=0
while cnt<3000:
    n=random.randint(4,12); p=random.uniform(0.3,0.9)
    e=[(i,j) for i,j in itertools.combinations(range(n),2) if random.random()<p]
    deg=[0]*n
    for i,j in e: deg[i]+=1; deg[j]+=1
    if min(deg)<2: continue
    test(n,e,'rand'); cnt+=1
print('tested',tested,'bad=',bad)
sys.exit(1 if bad else 0)
