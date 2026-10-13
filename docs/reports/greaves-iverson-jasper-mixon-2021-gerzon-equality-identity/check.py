#!/usr/bin/env python3
"""Gerzon equality in finite orthogonal geometries: does a^2 = (d+2) b hold without congruence conditions?

Source: Greaves, Iverson, Jasper, Mixon, "Frames over finite fields: Equiangular lines in orthogonal
geometry", arXiv:2012.13642v2, Example 5.8: a^2 = (d+2) b holds whenever Gerzon's bound is saturated,
d is not 0 mod p and d is not 1 mod p; "The authors do not know if these conditions can be removed."

For small (d, p) with d = 0 or 1 mod p (and a few controls with d not 0 or 1 mod p), this program enumerates every (a, b)-equiangular system of
N = d(d+1)/2 vectors with a^2 != b in F_p^d, for both isometry classes of nondegenerate symmetric
forms (diag(1,...,1) and diag(1,...,1,n) with n a nonsquare). It lists the systems up to sign and
reports whether a^2 = (d+2) b holds for each.  Exhaustive clique search; standard library only.
Prints one line per case and ALL_OK when every found system satisfies the identity.
"""
import itertools, sys

def nonsquare(p):
    sq = {x * x % p for x in range(1, p)}
    return next(x for x in range(2, p) if x not in sq)

def cliques_of_size(adj, k):
    # Bron-Kerbosch with pivoting, collecting maximal cliques of size >= k; returns cliques of size exactly k
    out = []
    def bk(R, P, X):
        if len(R) + len(P) < k: return
        if not P and not X:
            if len(R) >= k: out.append(R)
            return
        u = max(P | X, key=lambda w: len(adj[w] & P))
        for v in list(P - adj[u]):
            bk(R | {v}, P & adj[v], X & adj[v]); P = P - {v}; X = X | {v}
    bk(frozenset(), frozenset(adj), frozenset())
    res = set()
    for c in out:
        for s in itertools.combinations(sorted(c), k): res.add(s)
    return res

def run(d, p, diag):
    N = d * (d + 1) // 2
    B = lambda u, v: sum(diag[i] * u[i] * v[i] for i in range(d)) % p
    vecs = [v for v in itertools.product(range(p), repeat=d) if any(v)]
    found = []
    for a in range(p):
        lines = []
        seen = set()
        for v in vecs:
            if B(v, v) != a: continue
            m = tuple((-x) % p for x in v)
            if m in seen: continue
            seen.add(v); lines.append(v)
        for b in range(p):
            if (a * a - b) % p == 0: continue
            adj = {i: set() for i in range(len(lines))}
            for i, j in itertools.combinations(range(len(lines)), 2):
                if B(lines[i], lines[j]) ** 2 % p == b:
                    adj[i].add(j); adj[j].add(i)
            adj = {i: frozenset(s) for i, s in adj.items()}
            for c in cliques_of_size(adj, N):
                found.append((a, b, (a * a - (d + 2) * b) % p == 0))
    return found

ok = True
# exceptional classes d = 0 or 1 mod p first, then non-exceptional controls
for d, p in ((3, 3), (4, 3), (5, 5), (2, 3), (2, 5), (3, 5), (2, 7), (3, 7)):
    for label, diag in (("plus", [1] * d), ("minus", [1] * (d - 1) + [nonsquare(p)])):
        res = run(d, p, diag)
        pairs = sorted({(a, b) for a, b, _ in res})
        holds = all(h for _, _, h in res)
        ok &= holds
        print(f"d={d} p={p} (d mod p = {d % p}) form={label}: systems={len(res)} parameters={pairs} identity_holds={holds}")
        sys.stdout.flush()
print("ALL_OK" if ok else "FAIL")
sys.exit(0 if ok else 1)
