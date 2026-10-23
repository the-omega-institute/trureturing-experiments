#!/usr/bin/env python3
"""Checks for Q. Zhang and Y. Tian, "On the Steiner k-diameter and Steiner (k,k')-radius of trees", arXiv:2511.22492v1,
Conjecture 5.1: for a tree T of order n, n >= k >= 3, k > k' >= 1, Sd_k(T) <= k/(k-k') Sr_{k,k'}(T) - k'(k'-1)/(k-k'),
i.e. (k-k') Sd_k + k'(k'-1) <= k Sr_{k,k'}. Steiner distance d_T(S) = edges of the minimal subtree containing S
(an edge counts iff S meets both sides of it); Sd_k = max_{|S|=k} d_T(S); eps_{k,k'}(B) = max_{B <= S, |S|=k} d_T(S);
Sr_{k,k'} = min_{|B|=k'} eps_{k,k'}(B).
Fidelity: the equality cases of Theorems 4.4 and 4.6 (P_2(a,b;x): Sd_k = kx+1, Sr_{k,2} = (k-2)x+1; P_3(a,b;x):
Sd_k = kx+2, Sr_{k,3} = (k-3)x+2). Claim: every nonisomorphic tree with n <= 12 (counts A000055), every k, k';
also the pointwise strengthening k eps_{k,k'}(B) >= (k-k') Sd_k + k' d_T(B) for every k'-set B, and Theorem 3.4.
Exits 1 on failure.
"""
import itertools, sys

def trees(nmax):
    # nonisomorphic free trees via leaf extension + canonical form (AHU at centers)
    def canon(adj):
        n = len(adj)
        deg = [len(a) for a in adj]; leaves = [v for v in range(n) if deg[v] <= 1]; rem = n; layer = leaves
        removed = [False]*n
        while rem > 2:
            nl = []
            for v in layer:
                removed[v] = True; rem -= 1
                for u in adj[v]:
                    if not removed[u]:
                        deg[u] -= 1
                        if deg[u] == 1: nl.append(u)
            layer = nl
        centers = [v for v in range(n) if not removed[v]]
        def enc(v, p):
            return '(' + ''.join(sorted(enc(u, v) for u in adj[v] if u != p)) + ')'
        return min(enc(c, -1) for c in centers) if len(centers) == 1 else min(enc(centers[0], -1), enc(centers[1], -1), ''.join(sorted([enc(centers[0], centers[1]), enc(centers[1], centers[0])])))
    level = {canon([[]]): [[]]}
    out = {1: [[]]}
    cur = [[[]]]
    for n in range(2, nmax + 1):
        nxt = {}
        for adj in cur:
            for v in range(len(adj)):
                a = [list(x) for x in adj] + [[v]]; a[v].append(len(adj))
                c = canon(a)
                if c not in nxt: nxt[c] = a
        cur = list(nxt.values()); out[n] = cur
    return out
def sides(adj):
    n = len(adj); masks = []
    for v in range(n):
        for u in adj[v]:
            if u > v:
                # side of v when edge (v,u) removed
                seen = {v}; st = [v]
                while st:
                    x = st.pop()
                    for y in adj[x]:
                        if (x, y) in ((v, u), (u, v)): continue
                        if y not in seen: seen.add(y); st.append(y)
                masks.append(sum(1 << x for x in seen))
    return masks
def check_tree(adj, stats):
    n = len(adj); full = (1 << n) - 1; sm = sides(adj)
    d = [0]*(1 << n)
    for S in range(1, 1 << n):
        d[S] = sum(1 for m in sm if (S & m) and (S & ~m & full))
    bysize = [[] for _ in range(n + 1)]
    for S in range(1, 1 << n): bysize[bin(S).count('1')].append(S)
    for k in range(3, n + 1):
        D = max(d[S] for S in bysize[k])
        for r in range(1, k):
            RB = {}
            for S in bysize[k]:
                bits = [i for i in range(n) if S >> i & 1]
                for c in itertools.combinations(bits, r):
                    B = sum(1 << i for i in c)
                    if d[S] > RB.get(B, -1): RB[B] = d[S]
            R = min(RB.values())
            stats['conj'] += 1
            if (k - r) * D + r * (r - 1) > k * R: stats['conj_fail'].append((n, k, r))
            if (k - r) * D > k * R: stats['t34_fail'].append((n, k, r))
            for B, v in RB.items():
                stats['pw'] += 1
                if k * v < (k - r) * D + r * d[B]: stats['pw_fail'].append((n, k, r))
            if (k - r) * D + r * (r - 1) == k * R: stats['tight'][r] = stats['tight'].get(r, 0) + 1
def caterpillar(l, a, b, x):
    adj = [[] for _ in range(l)]
    for i in range(l - 1): adj[i].append(i + 1); adj[i + 1].append(i)
    for end, cnt in ((0, a), (l - 1, b)):
        for _ in range(cnt):
            prev = end
            for _ in range(x):
                v = len(adj); adj.append([prev]); adj[prev].append(v); prev = v
    return adj
def params(adj, k, r):
    n = len(adj); full = (1 << n) - 1; sm = sides(adj)
    dd = lambda S: sum(1 for m in sm if (S & m) and (S & ~m & full))
    D = max(dd(sum(1 << i for i in c)) for c in itertools.combinations(range(n), k))
    R = min(max(dd(sum(1 << i for i in set(B) | set(c))) for c in itertools.combinations([i for i in range(n) if i not in B], k - r)) for B in itertools.combinations(range(n), r))
    return D, R
ok = True
for (l, a, b, x, k) in ((2, 2, 2, 2, 4), (2, 3, 2, 1, 5), (2, 2, 2, 1, 4), (3, 2, 2, 2, 4), (3, 3, 3, 1, 6), (3, 3, 2, 1, 5)):
    D, R = params(caterpillar(l, a, b, x), k, l)
    want = (k * x + l - 1, (k - l) * x + l - 1)
    good = (D, R) == want and (k - l) * D + l * (l - 1) == k * R
    print(f"fidelity P_{l}({a},{b};{x}), k = {k}: Sd_k = {D}, Sr_k,{l} = {R}, source {want}, equality in Conjecture 5.1: {good}")
    ok &= good
N = 12
T = trees(N)
counts = [len(T[n]) for n in range(1, N + 1)]
print(f"nonisomorphic trees n = 1..{N}: {counts}")
ok &= counts == [1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551]
stats = {'conj': 0, 'conj_fail': [], 't34_fail': [], 'pw': 0, 'pw_fail': [], 'tight': {}}
for n in range(3, N + 1):
    for adj in T[n]:
        check_tree(adj, stats)
print(f"(T, k, k') instances: {stats['conj']}; Conjecture 5.1 failures: {len(stats['conj_fail'])}; Theorem 3.4 failures: {len(stats['t34_fail'])}; pointwise checks: {stats['pw']}, failures: {len(stats['pw_fail'])}; equality instances by k': {dict(sorted(stats['tight'].items()))}")
ok &= not stats['conj_fail'] and not stats['t34_fail'] and not stats['pw_fail']
print("ALL_CHECKS_OK" if ok else "CHECK_FAIL")
sys.exit(0 if ok else 1)
