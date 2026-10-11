#!/usr/bin/env python3
"""Ferenczi-Zamboni, arXiv:2507.17370v2, Section 4, Conjecture: "A language on k letters satisfying a symmetric
order condition and with no connection has palindromic complexity 1 for n even, k for n odd."
Languages are factorial and extendable (Section 1). Order condition (Def. 2.x): axb, cxd in L, a != c, b != d
=> (a <_A c iff b <_D d); symmetric: <_A is the reverse of <_D. Connection at a bispecial w: consecutive
a <_A a', consecutive b <_D b' with awb, a'wb' in L and neither awb' nor a'wb in L.

Family L_k (k >= 4): letters 1..k, c = k-1, d = k; edges (i,d) for 1<=i<=k-2, (c,j) for all j, (d,1);
L_k = all finite directed paths (and the empty word). <_D numerical, <_A reverse.
Checks for k = 4..8 and all words up to length 12 (exhaustive):
 (1) factorial, extendable on both sides, every letter used;
 (2) the order condition for every x of length <= 10 (not only bispecial ones);
 (3) the only bispecial word is the empty word, and it has no connection;
 (4) factor complexity p(n) = (k-1)n + 1;
 (5) palindromic complexity: PC(0)=1, PC(1)=k, PC(2m)=1, PC(2m+1)=3 < k for m >= 1."""
import sys, itertools
bad=0
def rep(l,ok):
    global bad; print(l,ok); bad+=0 if ok else 1
for k in range(4,9):
    c,d=k-1,k
    E={(i,d) for i in range(1,k-1)}|{(c,j) for j in range(1,k+1)}|{(d,1)}
    succ={a:[b for b in range(1,k+1) if (a,b) in E] for a in range(1,k+1)}
    pred={b:[a for a in range(1,k+1) if (a,b) in E] for b in range(1,k+1)}
    NMAX=12
    L={0:[()]}
    L[1]=[(a,) for a in range(1,k+1)]
    for n in range(2,NMAX+1):
        L[n]=[w+(b,) for w in L[n-1] for b in succ[w[-1]]]
    S={w for n in L for w in L[n]}
    ok=all(w[1:] in S and w[:-1] in S for n in range(1,NMAX+1) for w in L[n])
    ok&=all(succ[a] and pred[a] for a in range(1,k+1))
    rep(f'k={k} (1) factorial, extendable, all letters',ok)
    ok=True
    for n in range(0,NMAX-1):
        for x in L[n]:
            ext=[(a,b) for a in range(1,k+1) for b in range(1,k+1) if (a,)+x+(b,) in S]
            for (a,b),(cc,dd) in itertools.product(ext,ext):
                if a!=cc and b!=dd and ((a>cc)!=(b<dd)): ok=False
    rep(f'k={k} (2) symmetric order condition (|x|<=10)',ok)
    def left(x): return {a for a in range(1,k+1) if (a,)+x in S}
    def right(x): return {b for b in range(1,k+1) if x+(b,) in S}
    bis=[x for n in range(0,NMAX) for x in L[n] if len(left(x))>1 and len(right(x))>1]
    ok=bis==[()]
    # connection at the empty word: consecutive in <_A (reverse numerical): (a, a') with a' = a-1 means a <_A a'
    conn=False
    for a in range(2,k+1):
        ap=a-1  # a <_A ap, consecutive
        for b in range(1,k):
            bp=b+1  # b <_D bp, consecutive
            if (a,b) in E and (ap,bp) in E and (a,bp) not in E and (ap,b) not in E: conn=True
    ok&=not conn
    rep(f'k={k} (3) only bispecial is empty, no connection',ok)
    ok=all(len(L[n])==(k-1)*n+1 for n in range(0,NMAX+1))
    rep(f'k={k} (4) p(n)=(k-1)n+1',ok)
    pc=[sum(1 for w in L[n] if w==w[::-1]) for n in range(0,NMAX+1)]
    ok=pc[0]==1 and pc[1]==k and all(pc[n]==1 for n in range(2,NMAX+1,2)) and all(pc[n]==3 for n in range(3,NMAX+1,2))
    rep(f'k={k} (5) PC={pc}',ok)
print('bad=',bad); sys.exit(1 if bad else 0)
