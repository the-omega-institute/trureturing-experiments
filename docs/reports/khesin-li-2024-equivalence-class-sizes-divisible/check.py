#!/usr/bin/env python3
"""Khesin-Li (arXiv:2406.12083v1), Conjecture 7.1: for positive integers n > k, the equivalence classes of [n,k] codes
with a single connected component have sizes divisible by 3^k.

(1) Reproduces the Section 7 tabulation (normalized encoders; one operation of Conjecture 6.3 followed by renormalization:
    drop input-input edges, RREF of the input-output matrix, pivots moved first; results with pivot-pivot edges are not
    merged; disjoint sets) and compares the connected class sizes with the paper's Figures for [2,1], [3,1], [4,1], [5,1],
    [3,2], [4,2].
(2) Checks the proof's group action T_i = LC(p_i) o LC(u_i) on every normalized encoder for (n,k) up to (5,3): each local
    complementation stays normalized, T_i^3 = id, T_i T_j = T_j T_i, connectivity is preserved, and no nonidentity element
    of (Z/3)^k fixes a connected normalized encoder.
Exit 0 and RESULT PASS when every check holds.
"""
import itertools, sys
from collections import defaultdict
def run(n,k,pp_rule):
    m=n-k; N=k+n
    # vertices: 0..k-1 inputs, k..k+k-1 pivots, k+k.. free
    def normal_graphs():
        free_pairs=[(a,b) for a in range(m) for b in range(a+1,m)]
        for B in itertools.product([0,1],repeat=k*m):
            for C in itertools.product([0,1],repeat=k*m):
                for D in itertools.product([0,1],repeat=len(free_pairs)):
                    A=[[0]*N for _ in range(N)]
                    for i in range(k): A[i][k+i]=A[k+i][i]=1
                    for i in range(k):
                        for a in range(m):
                            if B[i*m+a]: A[i][2*k+a]=A[2*k+a][i]=1
                            if C[i*m+a]: A[k+i][2*k+a]=A[2*k+a][k+i]=1
                    for t,(a,b) in enumerate(free_pairs):
                        if D[t]: A[2*k+a][2*k+b]=A[2*k+b][2*k+a]=1
                    yield A
    key=lambda A: tuple(tuple(r) for r in A)
    def lc(A,v):
        A=[r[:] for r in A]; nb=[u for u in range(N) if A[v][u]]
        for x,y in itertools.combinations(nb,2): A[x][y]^=1; A[y][x]^=1
        return A
    def normalize(A):
        A=[r[:] for r in A]
        for i in range(k):
            for j in range(k):
                if i!=j: A[i][j]=0
        # RREF of k x n input-output matrix (columns = outputs k..N-1)
        M=[[A[i][k+c] for c in range(n)] for i in range(k)]
        piv=[]; r=0
        for c in range(n):
            pr=next((i for i in range(r,k) if M[i][c]),None)
            if pr is None: continue
            M[r],M[pr]=M[pr],M[r]
            for i in range(k):
                if i!=r and M[i][c]: M[i]=[x^y for x,y in zip(M[i],M[r])]
            piv.append(c); r+=1
            if r==k: break
        if r<k: return None
        for i in range(k):
            for c in range(n): A[i][k+c]=A[k+c][i]=M[i][c]
        # output permutation: pivots first in order, then the rest in order
        order=piv+[c for c in range(n) if c not in piv]
        perm=list(range(k))+[k+c for c in order]
        B=[[A[perm[x]][perm[y]] for y in range(N)] for x in range(N)]
        if any(B[k+i][k+j] for i in range(k) for j in range(k) if i!=j):
            if pp_rule=='skip': return None
        return B
    G=list(normal_graphs()); idx={key(A):t for t,A in enumerate(G)}
    par=list(range(len(G)))
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    def union(a,b): par[find(a)]=find(b)
    for t,A in enumerate(G):
        res=[]
        for v in range(N): res.append(normalize(lc(A,v)))
        for a,b in itertools.combinations(range(m),2):   # free output transpositions
            perm=list(range(N)); perm[2*k+a],perm[2*k+b]=perm[2*k+b],perm[2*k+a]
            res.append(normalize([[A[perm[x]][perm[y]] for y in range(N)] for x in range(N)]))
        for i,j in itertools.permutations(range(k),2):    # row operation
            Bm=[r[:] for r in A]
            for c in range(n): Bm[i][k+c]^=A[j][k+c]; Bm[k+c][i]=Bm[i][k+c]
            res.append(normalize(Bm))
        for R in res:
            if R is not None and key(R) in idx: union(t,idx[key(R)])
    def connected(A):
        seen={0}; st=[0]
        while st:
            x=st.pop()
            for y in range(N):
                if A[x][y] and y not in seen: seen.add(y); st.append(y)
        return len(seen)==N
    cls=defaultdict(list)
    for t in range(len(G)): cls[find(t)].append(t)
    sizes=sorted(len(v) for v in cls.values() if all(connected(G[t]) for t in v))
    mixed=sum(1 for v in cls.values() if len(set(connected(G[t]) for t in v))>1)
    return sizes,mixed
def check(n,k):
    m=n-k; N=k+n; fp=[(a,b) for a in range(m) for b in range(a+1,m)]
    def graphs():
        for B in itertools.product([0,1],repeat=k*m):
            for C in itertools.product([0,1],repeat=k*m):
                for D in itertools.product([0,1],repeat=len(fp)):
                    A=[[0]*N for _ in range(N)]
                    for i in range(k): A[i][k+i]=A[k+i][i]=1
                    for i in range(k):
                        for a in range(m):
                            if B[i*m+a]: A[i][2*k+a]=A[2*k+a][i]=1
                            if C[i*m+a]: A[k+i][2*k+a]=A[2*k+a][k+i]=1
                    for t,(a,b) in enumerate(fp):
                        if D[t]: A[2*k+a][2*k+b]=A[2*k+b][2*k+a]=1
                    yield A
    key=lambda A: tuple(map(tuple,A))
    def lc(A,v):
        A=[r[:] for r in A]; nb=[u for u in range(N) if A[v][u]]
        for x,y in itertools.combinations(nb,2): A[x][y]^=1; A[y][x]^=1
        return A
    def normalized(A):
        if any(A[i][j] for i in range(k) for j in range(k)): return False
        if any(A[i][k+j]!=(1 if i==j else 0) for i in range(k) for j in range(k)): return False
        if any(A[k+i][k+j] for i in range(k) for j in range(k)): return False
        return all(A[x][x]==0 for x in range(N))
    def connected(A):
        seen={0}; st=[0]
        while st:
            x=st.pop()
            for y in range(N):
                if A[x][y] and y not in seen: seen.add(y); st.append(y)
        return len(seen)==N
    T=lambda A,i: lc(lc(A,i),k+i)
    bad=0; cnt=0
    for A in graphs():
        conn=connected(A)
        for i in range(k):
            Ti=T(A,i)
            if not normalized(lc(A,i)) or not normalized(Ti): bad+=1
            if key(T(T(Ti,i),i))!=key(A): bad+=1
            for j in range(k):
                if key(T(T(A,i),j))!=key(T(T(A,j),i)): bad+=1
        if conn:
            cnt+=1
            for e in itertools.product(range(3),repeat=k):
                if any(e):
                    X=A
                    for i in range(k):
                        for _ in range(e[i]): X=T(X,i)
                    if key(X)==key(A): bad+=1
            if any(connected(T(A,i))!=conn for i in range(k)): bad+=1
    return cnt,bad
PAPER = {(2, 1): [3], (3, 1): [3, 21], (4, 1): [3, 30, 45, 54, 84, 198],
         (5, 1): [3, 39, 78, 84, 84, 204, 297, 306, 315, 360, 540, 558, 1332, 1404, 2376, 3024, 3276],
         (3, 2): [9], (4, 2): [36, 45, 99, 234]}
ok = True
for nk, want in PAPER.items():
    s, mx = run(*nk, 'skip')
    good = s == sorted(want) and mx == 0 and all(x % 3 ** nk[1] == 0 for x in s)
    print(nk, 'class sizes', s, 'match paper:', s == sorted(want), 'mixed classes', mx); ok = ok and good
for nk in [(2, 1), (3, 1), (4, 1), (5, 1), (3, 2), (4, 2), (5, 2), (4, 3), (5, 3)]:
    cnt, bad = check(*nk); print(nk, 'connected normalized encoders', cnt, 'action failures', bad); ok = ok and bad == 0
print('RESULT PASS' if ok else 'RESULT FAIL'); sys.exit(0 if ok else 1)
