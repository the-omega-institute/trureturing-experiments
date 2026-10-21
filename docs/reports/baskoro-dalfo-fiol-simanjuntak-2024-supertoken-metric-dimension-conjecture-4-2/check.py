"""Checks for Conjecture 4.2 of Baskoro, Dalfó, Fiol and Simanjuntak, "On some metric properties of supertoken graphs",
arXiv:2412.20558v1: dim(F_k(K_n)) = n - 1. Vertices: x in N^n with sum k; adjacent when one token moves between two
distinct vertices of K_n. Distances by breadth-first search on the adjacency (not by a formula).
(a) Reproduce Proposition 4.1 (the n-1 resolving set) and the n = 3, 4 equality of Proposition 4.3 for small k.
(b) The landmark set W = (12, 13, 14, 56, 57, 58) resolves F_2(K_8): dim(F_2(K_8)) <= 6 < 7.
(c) The family: the 3r star landmarks resolve F_2(K_{4r}) for r = 2, 3, 4.
(d) Exact dim(F_2(K_n)) for small n by exhaustive search."""
import itertools, sys
from collections import deque

def vertices(n, k):
    out = []
    def rec(i, left, cur):
        if i == n - 1:
            out.append(tuple(cur + [left])); return
        for a in range(left + 1):
            rec(i + 1, left - a, cur + [a])
    rec(0, k, [])
    return out

def neighbours(x):
    n = len(x)
    for i in range(n):
        if x[i] == 0: continue
        for j in range(n):
            if j == i: continue
            y = list(x); y[i] -= 1; y[j] += 1
            yield tuple(y)

def all_distances(V):
    idx = {v: t for t, v in enumerate(V)}
    D = []
    for s in V:
        d = [-1] * len(V); d[idx[s]] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for w in neighbours(u):
                if d[idx[w]] < 0:
                    d[idx[w]] = d[idx[u]] + 1; q.append(w)
        D.append(d)
    return idx, D

def resolves(V, D, W):
    seen = set()
    for t in range(len(V)):
        key = tuple(D[w][t] for w in W)
        if key in seen: return False
        seen.add(key)
    return True

def metric_dimension(V, D, upto):
    for s in range(1, upto + 1):
        for W in itertools.combinations(range(len(V)), s):
            if resolves(V, D, W): return s
    return None

def pair(n, i, j):
    x = [0] * n; x[i - 1] += 1; x[j - 1] += 1; return tuple(x)

def main():
    bad = 0
    print("(a) Proposition 4.1 and Proposition 4.3 (n = 3, 4)")
    for n, k in [(3, 1), (3, 2), (3, 3), (3, 4), (4, 1), (4, 2), (4, 3)]:
        V = vertices(n, k); idx, D = all_distances(V)
        C = [idx[tuple(k if t == i else 0 for t in range(n))] for i in range(n - 1)]
        prop41 = resolves(V, D, C)
        dim = metric_dimension(V, D, n - 1)
        ok = prop41 and dim == n - 1
        bad += not ok
        print(f"  n={n} k={k}: |V|={len(V)} Prop 4.1 set resolves {prop41}; dim={dim} (conjectured {n - 1}) {'OK' if ok else 'MISMATCH'}")
    print("(b) F_2(K_8) with W = (12,13,14,56,57,58)")
    n = 8; V = vertices(n, 2); idx, D = all_distances(V)
    W = [idx[pair(n, *p)] for p in [(1, 2), (1, 3), (1, 4), (5, 6), (5, 7), (5, 8)]]
    r8 = resolves(V, D, W); bad += not r8
    print(f"  |V|={len(V)} edges={sum(len(list(neighbours(v))) for v in V)//2} diameter={max(max(r) for r in D)} resolves={r8} -> dim <= 6 < 7")
    print("(c) the star family on F_2(K_{4r})")
    for r in [2, 3, 4]:
        n = 4 * r; V = vertices(n, 2); idx, D = all_distances(V)
        W = [idx[pair(n, 4 * b - 3, 4 * b - c)] for b in range(1, r + 1) for c in (2, 1, 0)]
        ok = resolves(V, D, W); bad += not ok
        print(f"  r={r} n={n}: {len(W)} landmarks resolve {ok}; conjectured n-1={n - 1}")
    print("(d) exact dim(F_2(K_n)), small n")
    for n in range(3, 8):
        V = vertices(n, 2); idx, D = all_distances(V)
        dim = metric_dimension(V, D, n - 1)
        print(f"  n={n}: dim={dim} (conjectured {n - 1}){'  <- counterexample' if dim is not None and dim < n - 1 else ''}")
    print("(e) exact dim(F_2(K_8))")
    V = vertices(8, 2); idx, D = all_distances(V)
    dim8 = metric_dimension(V, D, 7)
    bad += dim8 != 6
    print(f"  n=8: dim={dim8} (conjectured 7)")
    print(f"EXIT={1 if bad else 0}")
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
