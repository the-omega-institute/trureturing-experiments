"""arXiv:2307.01648 closing conjecture: Algorithm 1 (proof of Lemma 4) finds a shortest path in G(P).

Algorithm 1: relabel so that P[1] >= P[2] >= ...; for each x != 1 in label order, pair the k-th occurrence of x in
the current word with the k-th occurrence of x in v, then apply all those swaps. Every admissible run (every
reference symbol of maximal count, every order of the other symbols with non-increasing count) is tried; only
nontrivial swaps (i != j) are counted. Distances are breadth-first search distances in the 2-swap graph.
Checks:
  1. family F_k: w = 1^k 2^k, v = 2 1^(k-1) 2^(k-1) 1, for 2 <= k <= 7: distance 1 and every admissible run
     generates k nontrivial swaps. The stage swaps are applied in the order of k and, separately, in reverse order;
     for each run at least one of the two orders is a path of legal 2-swaps from w to v;
  2. family G_m (search seat): w = 1 2^m 1^m, v = 2^m 1^(m+1), for 2 <= m <= 6: distance 1, Algorithm 1 uses m swaps;
  3. exhaustive binary census n <= 9 and ternary census n <= 6: the number of pairs where no admissible run is
     shortest, and the smallest n with such a pair.
"""
from collections import deque
from itertools import permutations, product

def occ(w, x):
    return [i for i, c in enumerate(w) if c == x]

def run(w, v, order, reverse=False):
    cur = list(w); count = 0; legal = True
    for x in order:
        S = list(zip(occ(cur, x), occ(v, x)))
        for i, j in (S[::-1] if reverse else S):
            if i != j:
                if cur[i] == cur[j]:
                    legal = False
                cur[i], cur[j] = cur[j], cur[i]; count += 1
    return count, tuple(cur) == tuple(v), legal

def admissible_orders(w):
    syms = sorted(set(w)); cnt = {s: w.count(s) for s in syms}; mx = max(cnt.values())
    for ref in [s for s in syms if cnt[s] == mx]:
        rest = [s for s in syms if s != ref]
        for o in permutations(rest):
            if all(cnt[o[t]] >= cnt[o[t + 1]] for t in range(len(o) - 1)):
                yield list(o)

def dist(w, v):
    w, v = tuple(w), tuple(v); d = {w: 0}; q = deque([w])
    while q:
        x = q.popleft()
        if x == v:
            return d[x]
        l = list(x)
        for i in range(len(l)):
            for j in range(i + 1, len(l)):
                if l[i] != l[j]:
                    u = l[:]; u[i], u[j] = u[j], u[i]; u = tuple(u)
                    if u not in d:
                        d[u] = d[x] + 1; q.append(u)

ok = True
for k in range(2, 8):
    w = [1] * k + [2] * k; v = [2] + [1] * (k - 1) + [2] * (k - 1) + [1]
    orders = list(admissible_orders(w))
    runs = [(o, run(w, v, o), run(w, v, o, reverse=True)) for o in orders]
    d = dist(w, v)
    good = d == 1 and len(runs) == 2 and all(f[0] == k and r[0] == k for _, f, r in runs) \
        and all(f[1:] == (True, True) or r[1:] == (True, True) for _, f, r in runs)
    ok &= good
    print(f"F_{k}: dist={d} " + "; ".join(f"order={o} forward(count,reaches,legal)={f} reverse={r}" for o, f, r in runs)
          + f" {'ok' if good else 'FAIL'}")
for m in range(2, 7):
    w = [1] + [2] * m + [1] * m; v = [2] * m + [1] * (m + 1)
    runs = [run(w, v, o) for o in admissible_orders(w)]
    d = dist(w, v)
    good = d == 1 and runs == [(m, True, True)]
    ok &= good
    print(f"G_{m}: dist={d} runs={runs} {'ok' if good else 'FAIL'}")
for alph, N in ((2, 9), (3, 6)):
    for n in range(2, N + 1):
        bad = total = 0
        for w in product(range(alph), repeat=n):
            for v in set(permutations(w)):
                total += 1
                d = dist(w, v)
                best = min(run(w, v, o)[0] for o in admissible_orders(list(w)))
                if best != d:
                    bad += 1
        print(f"alphabet {alph}, n={n}: pairs={total} pairs with no shortest admissible run={bad}")
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
