"""B. A. Rather, arXiv:2603.06930v3, Section 8 item 1: real-rootedness of the general position polynomial.

For every labelled simple graph on n <= N vertices the program computes
    psi(G) = sum_k alpha_k(G) x^k,  alpha_k = number of general position sets of size k,
where S is in general position if no three distinct vertices x, y, z of S have y on a shortest x-z path
(d(x,z) = d(x,y) + d(y,z) with all three distances finite), counts the real roots of psi(G) exactly
(sympy real_roots, with multiplicity) and checks
1. connected G: psi(G) has only real zeros iff G is complete;
2. all G: psi(G) has only real zeros iff every component of G is complete;
3. the Newton step: alpha_0 = 1, alpha_1 = n, alpha_2 = n(n-1)/2 for every G, and a real-rooted psi has degree n.
Prints RESULT PASS or RESULT FAIL.
"""
import itertools
import sys
from collections import deque

import sympy

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
x = sympy.symbols("x")


def distances(n, adj):
    D = [[None] * n for _ in range(n)]
    for s in range(n):
        D[s][s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if D[s][v] is None:
                    D[s][v] = D[s][u] + 1
                    q.append(v)
    return D


def psi(n, adj):
    D = distances(n, adj)
    c = [0] * (n + 1)
    for mask in range(1 << n):
        S = [i for i in range(n) if mask >> i & 1]
        gp = True
        for a, b, d in itertools.permutations(S, 3):
            if None not in (D[a][d], D[a][b], D[b][d]) and D[a][d] == D[a][b] + D[b][d]:
                gp = False
                break
        if gp:
            c[len(S)] += 1
    return c


def real_rooted(c):
    p = sympy.Poly(sum(ci * x ** i for i, ci in enumerate(c)), x)
    return len(p.real_roots()) == p.degree()


def connected(n, adj):
    seen, q = {0}, [0]
    while q:
        u = q.pop()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                q.append(v)
    return len(seen) == n


tot = conn = 0
fail = {"connected": 0, "cluster": 0, "low": 0, "degree": 0}
for n in range(1, N + 1):
    pairs = list(itertools.combinations(range(n), 2))
    for m in range(1 << len(pairs)):
        adj = [set() for _ in range(n)]
        E = [pairs[i] for i in range(len(pairs)) if m >> i & 1]
        for a, b in E:
            adj[a].add(b)
            adj[b].add(a)
        c = psi(n, adj)
        rr = real_rooted(c)
        tot += 1
        if c[0] != 1 or c[1] != n or (n >= 2 and c[2] != n * (n - 1) // 2):
            fail["low"] += 1
        if rr and max(k for k in range(n + 1) if c[k]) != n:
            fail["degree"] += 1
        cluster = all(b in adj[a] for v in range(n) for a in adj[v] for b in adj[v] if a != b)
        if rr != cluster:
            fail["cluster"] += 1
        if connected(n, adj):
            conn += 1
            if rr != (len(E) == len(pairs)):
                fail["connected"] += 1
ok = all(v == 0 for v in fail.values())
print("labelled graphs %d (connected %d) on 1..%d vertices; failures %s" % (tot, conn, N, fail))
print("RESULT " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
