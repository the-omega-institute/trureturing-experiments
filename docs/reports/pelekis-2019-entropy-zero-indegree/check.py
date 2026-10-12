"""Check of a proposed refutation of Pelekis, 'On the Shannon entropy of the number of vertices with zero in-degree in
randomly oriented hypergraphs' (arXiv:1808.04166v1; Le Matematiche 74(1) (2019) 119-130), Conjecture 1.3: an entropy
maximiser over D_{n,m,r} (r-uniform multi-hypergraphs, n vertices, m edges, parallel edges allowed) has all degrees
within one of each other.  X = number of vertices chosen by at least one edge (each edge picks a uniform endpoint).

Exhaustive over loopless multigraphs (r = 2): every multiset of m edges on n vertices, exact distribution of X by
enumerating all 2^m orientations, entropy in exact-rational-input floats; report the maximisers and whether any is
degree-balanced.  (n, m) = (3, 3..10), (4, 4..7), (5, 5), (5, 6), (6, 6).
"""
import itertools, math
from collections import Counter
bad = 0
def entropy(counts, total):
    return -sum(c/total*math.log2(c/total) for c in counts if c)
def dist(n, edges):
    m = len(edges); cnt = Counter()
    for choice in itertools.product((0, 1), repeat=m):
        sel = {e[c] for e, c in zip(edges, choice)}
        cnt[len(sel)] += 1
    return cnt
for (n, m) in [(3, k) for k in range(3, 11)] + [(4, k) for k in range(4, 8)] + [(5, 5), (5, 6), (6, 6)]:
    pairs = list(itertools.combinations(range(n), 2))
    best = -1; maxis = []
    for ms in itertools.combinations_with_replacement(range(len(pairs)), m):
        edges = [pairs[i] for i in ms]
        H = entropy(dist(n, edges).values(), 2 ** m)
        if H > best + 1e-12: best = H; maxis = [edges]
        elif abs(H - best) <= 1e-12: maxis.append(edges)
    def degs(edges):
        d = [0] * n
        for a, b in edges: d[a] += 1; d[b] += 1
        return d
    balanced = [e for e in maxis if max(degs(e)) - min(degs(e)) <= 1]
    print(f"n={n} m={m} maxH={best:.6f} maximisers={len(maxis)} balanced_maximisers={len(balanced)} example_degrees={sorted(degs(maxis[0]), reverse=True)}")
    if balanced: bad += 1  # the refutation predicts no balanced maximiser in these cases
print("bad=", bad)
