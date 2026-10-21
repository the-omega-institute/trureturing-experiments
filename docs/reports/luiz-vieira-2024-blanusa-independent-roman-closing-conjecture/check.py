"""Luiz-Vieira, arXiv:2406.11688v1 (DMGT 45 (2025)), Section 6: gamma_[kR](G) = i_[kR](G) for k >= 2 fails on the
generalized Blanusa snarks.

The script builds B_i^t (t = 1, 2) from the blocks B_0^1, B_0^2, L read off Figure 5 and the recursion of Section 4
(B_1^t, B_2^t as listed; B_i^t from B_{i-2}^t and the link graph LG_i). It checks that each B_i^t is cubic of order
8i + 10 and has girth 5, that B_1^t and B_2^t are not 3-edge-colourable, and that the domination and independent
domination numbers equal the source's Theorem 19 (2i + 4 for t = 1 and odd i >= 3, else 2i + 3) for i <= 3. For
k = 8i + 13 the uniform labelling r = 2i + 4 is a [k]-Roman dominating function of weight (8i+10)(2i+4), while the
minimum weight of an independent [k]-RDF is (k+1)(2i+3) = (8i+10)(2i+4) + 2 (exhaustive for i <= 2).
"""
import collections
import itertools
import sys

B1 = [(0, 1), (1, 2), (2, 3), (3, 9), (9, 8), (8, 7), (7, 6), (6, 4), (4, 0), (1, 7), (2, 5), (5, 8), (4, 5)]
B2 = [(0, 1), (1, 3), (1, 2), (2, 8), (2, 5), (0, 4), (4, 6), (4, 5), (5, 7), (3, 7), (7, 9), (6, 8), (8, 9)]
LE = [(0, 1), (1, 2), (2, 4), (4, 7), (7, 6), (6, 5), (5, 3), (3, 0), (1, 6), (3, 4)]
a, b, c, d = ('B', 0), ('B', 6), ('B', 3), ('B', 9)
X = lambda j: ('L', j, 0)
W = lambda j: ('L', j, 2)
Y = lambda j: ('L', j, 5)
Z = lambda j: ('L', j, 7)
E2 = lambda u, v: frozenset((u, v))


def block(t):
    return {E2(('B', u), ('B', v)) for u, v in (B1 if t == 1 else B2)}


def link(j):
    return {E2(('L', j, u), ('L', j, v)) for u, v in LE}


def edges(t, i):
    if i == 1:
        return block(t) | link(1) | {E2(c, Y(1)), E2(d, X(1)), E2(a, Z(1)), E2(b, W(1))}
    if i == 2:
        return block(t) | link(1) | link(2) | {E2(c, Y(1)), E2(d, X(1)), E2(W(1), Y(2)), E2(Z(1), X(2)),
                                                E2(a, Z(2)), E2(b, W(2))}
    prev = edges(t, i - 2) - {E2(a, Z(i - 2)), E2(b, W(i - 2))}
    lg = link(i - 1) | link(i) | {E2(W(i - 1), Y(i)), E2(Z(i - 1), X(i))}
    return prev | lg | {E2(W(i - 2), Y(i - 1)), E2(Z(i - 2), X(i - 1)), E2(a, Z(i)), E2(b, W(i))}


def graph(t, i):
    E = edges(t, i)
    adj = collections.defaultdict(set)
    for e in E:
        u, v = tuple(e)
        adj[u].add(v); adj[v].add(u)
    return adj


def girth(adj):
    g = 99
    for s in adj:
        dist, par, q = {s: 0}, {s: None}, collections.deque([s])
        while q:
            u = q.popleft()
            for w in adj[u]:
                if w not in dist:
                    dist[w] = dist[u] + 1; par[w] = u; q.append(w)
                elif par[u] != w:
                    g = min(g, dist[u] + dist[w] + 1)
    return g


def three_edge_colourable(adj):
    E = sorted({tuple(sorted(map(str, (u, v)))) + (u, v) for u in adj for v in adj[u] if str(u) < str(v)})
    E = [(u, v) for *_, u, v in E]
    col = {}
    used = collections.defaultdict(set)

    def bt(k):
        if k == len(E):
            return True
        u, v = E[k]
        for cc in range(3):
            if cc not in used[u] and cc not in used[v]:
                used[u].add(cc); used[v].add(cc)
                if bt(k + 1):
                    return True
                used[u].discard(cc); used[v].discard(cc)
        return False
    return bt(0)


def min_dominating(adj, independent):
    V = list(adj)
    for m in range(1, len(V) + 1):
        for S in itertools.combinations(V, m):
            Ss = set(S)
            if independent and any(w in Ss for u in S for w in adj[u]):
                continue
            if all(v in Ss or adj[v] & Ss for v in V):
                return m
    return None


def min_indep_roman(adj, k):
    V = list(adj); n = len(V); best = None
    for m in range(1, n + 1):
        if best is not None and m * k > best:
            break
        for S in itertools.combinations(V, m):
            Ss = set(S)
            if any(w in Ss for u in S for w in adj[u]) or any(v not in Ss and not adj[v] & Ss for v in V):
                continue
            for lab in itertools.product((k, k + 1), repeat=m):
                L = dict(zip(S, lab))
                if all(sum(L.get(u, 0) for u in adj[v]) >= k + len(adj[v] & Ss) for v in V if v not in Ss):
                    wgt = sum(lab)
                    best = wgt if best is None or wgt < best else best
    return best


ok = True
for t in (1, 2):
    for i in range(1, 7):
        adj = graph(t, i)
        n = len(adj)
        cubic = all(len(adj[v]) == 3 for v in adj)
        good = cubic and n == 8 * i + 10 and girth(adj) == 5
        line = f"B_{i}^{t}: order {n} cubic {cubic} girth {girth(adj)}"
        if i <= 2:
            col = three_edge_colourable(adj)
            good &= not col
            line += f" 3-edge-colourable {col}"
        if i <= 3:
            gam, ind = min_dominating(adj, False), min_dominating(adj, True)
            th19 = 2 * i + 4 if (t == 1 and i >= 3 and i % 2 == 1) else 2 * i + 3
            good &= gam == ind == th19
            line += f" gamma {gam} i {ind} (Theorem 19: {th19})"
        k = 8 * i + 13; r = 2 * i + 4
        unif = all(4 * r >= k + 3 for v in adj) and r < k
        line += f"; k = {k}: uniform r = {r} valid {unif}, weight {n * r}"
        good &= unif
        if i <= 2:
            ir = min_indep_roman(adj, k)
            good &= ir == (k + 1) * (2 * i + 3) == n * r + 2
            line += f", i_[kR] = {ir}"
        ok &= good
        print(line, "ok" if good else "FAIL")
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
