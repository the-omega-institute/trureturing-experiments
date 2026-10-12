"""Check of the Section 4 conjecture of Pelekis (arXiv:1808.04166v1; Le Matematiche 74(1) (2019) 119-130): "We
conjecture that circular hypergraphs are such that H(X_{C_{n,r}}) >= H(X_H), for all H in D_{n,n,r}" (r >= 3, n > r),
where C_{n,r} has as edges all n paths of size r of the cycle Z_n.  Exhaustive over r-uniform multi-hypergraphs with
n vertices and n edges; X = number of vertices with non-zero in-degree when each edge picks a uniform element.
For (n, r) = (4, 3), (5, 3), (6, 3), (5, 4) the maximum entropy exceeds H(X_{C_{n,r}}); an exact integer comparison
(sum of c log c over orientation counts) is printed for (4, 3).
"""
import itertools, math
from collections import Counter
bad = 0
def dist(edges, r):
    cnt = Counter()
    for ch in itertools.product(range(r), repeat=len(edges)):
        cnt[len({e[c] for e, c in zip(edges, ch)})] += 1
    return cnt
def H(cnt, tot): return -sum(c/tot*math.log2(c/tot) for c in cnt.values() if c)
for (n, r) in [(4, 3), (5, 3), (6, 3), (5, 4)]:
    sets = list(itertools.combinations(range(n), r)); tot = r ** n
    best, arg = -1.0, None
    for ms in itertools.combinations_with_replacement(range(len(sets)), n):
        e = [sets[i] for i in ms]; h = H(dist(e, r), tot)
        if h > best + 1e-12: best, arg = h, e
    circ = [tuple(sorted((i + j) % n for j in range(r))) for i in range(n)]
    dc = dist(circ, r); hc = H(dc, tot)
    print(f"n={n} r={r} max={best:.6f} circular={hc:.6f} argmax_edges={arg} counts_max={dict(dist(arg, r))} counts_circ={dict(dc)}")
    if not best > hc + 1e-9: bad += 1
    if (n, r) == (4, 3):
        # same total 81: H(A) > H(B)  iff  prod_B c^c > prod_A c^c  (exact integers)
        pa = 1
        for c in dist(arg, r).values(): pa *= c ** c
        pb = 1
        for c in dc.values(): pb *= c ** c
        print("exact (4,3): prod_circ c^c > prod_max c^c :", pb > pa)
        if not pb > pa: bad += 1
print("bad=", bad)
