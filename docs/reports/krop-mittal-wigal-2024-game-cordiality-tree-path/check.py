#!/usr/bin/env python3
"""Krop–Mittal–Wigal arXiv:2403.18060v1, Conjecture (Section 3): for any tree T of order n, c_g(T) <= c_g(P_n).
Exact minimax of the cordiality game (Admirable starts, labels 0, minimizes |e_1-e_0|; Impish labels 1, maximizes),
by memoized backward induction over (A, B) = (Admirable's vertices, Impish's vertices).
Checks: c_g(P_10) = 1; c_g(T) = 3 for the subdivided claw with arms 7, 1, 1 (edges (i,i+1), i<7, and (0,8), (0,9));
the Impish pairing strategy (0,4),(1,5),(2,6),(3,7),(8,9) forces |e_1-e_0| >= 3 against every Admirable play;
and the family T_r (2r leaves at vertex 0 of P_8) for r = 1, 2, 3 against the path of the same order."""
import sys
from functools import lru_cache
def cg(E,n):
    full=(1<<n)-1; m=len(E)
    @lru_cache(None)
    def F(A,B):
        free=full^(A|B)
        if free==0:
            q=sum(((A>>u)^(A>>v))&1 for u,v in E)
            return abs(2*q-m)
        a_turn=bin(A).count('1')==bin(B).count('1')
        vals=[F(A|(1<<v),B) if a_turn else F(A,B|(1<<v)) for v in range(n) if free&(1<<v)]
        return min(vals) if a_turn else max(vals)
    return F(0,0)
def pairing_lower(E,n,pairs):
    # Impish answers each Admirable move with its partner; minimum over all Admirable move orders of the final discrepancy
    mate={}
    for a,b in pairs: mate[a]=b; mate[b]=a
    full=(1<<n)-1; m=len(E)
    @lru_cache(None)
    def G(A,B):
        free=full^(A|B)
        if free==0:
            q=sum(((A>>u)^(A>>v))&1 for u,v in E); return abs(2*q-m)
        return min(G(A|(1<<v),B|(1<<mate[v])) for v in range(n) if free&(1<<v))
    return G(0,0)
bad=0
P10=[(i,i+1) for i in range(9)]; T=[(i,i+1) for i in range(7)]+[(0,8),(0,9)]
a,b=cg(P10,10),cg(T,10); print('cg(P10)=',a,'cg(T)=',b); bad+= (a!=1)+(b!=3)
lb=pairing_lower(T,10,[(0,4),(1,5),(2,6),(3,7),(8,9)]); print('pairing lower bound on T =',lb); bad+=(lb<3)
for r in (1,2):
    n=8+2*r; Tr=[(i,i+1) for i in range(7)]+[(0,8+k) for k in range(2*r)]; Pn=[(i,i+1) for i in range(n-1)]
    x,y=cg(Tr,n),cg(Pn,n); print(f'r={r} n={n} cg(T_r)={x} cg(P_n)={y}'); bad+=(x<=y)
print('bad=',bad); sys.exit(1 if bad else 0)
