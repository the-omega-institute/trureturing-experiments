"""Check of the conjecture after Eq. (87) in the Supplemental Material of Lee, Hhan, Cho and Kwon (arXiv:2507.17871v3):
chi'(H) <= O(t log t) with failure probability smaller than eps, for alpha i.i.d. uniform m-subsets of [k].

(a) Exact enumeration for small (k, m, alpha): chi'(H) by brute force, the lower bound Pr[chi' >= r] >= (m/k)^r of
Lemma 1, and the source's greedy bound (77) chi'(H) <= m (Delta(H) - 1) + 1.  (b) The refuting family t = 2^q,
eps = 2^-(t^2), c = 1, beta = 3: m = q, k = 3 t^2 + 6 q, alpha = 3 t^3 are integers, and for L >= 1 and q >= 32 (L + 1)
the event "the first r = L t q + 1 edges contain vertex 1" has probability (q/k)^r > eps, checked exactly as
(2q + 3) r < t^2 with q/k >= 2^-(2q+3), for L = 1..4 at the least admissible q.
"""
import itertools
from fractions import Fraction as F


def chromatic_index(edges):
    n = len(edges)
    for c in range(0, n + 1):
        def ok(i, col):
            if i == n:
                return True
            for x in range(c):
                if all(not (col[j] == x and (edges[j] & edges[i])) for j in range(i)):
                    col.append(x)
                    if ok(i + 1, col):
                        return True
                    col.pop()
            return False
        if ok(0, []):
            return c
    return n


okall = True
for k, m, alpha in [(4, 2, 3), (5, 2, 3), (5, 3, 3), (6, 2, 3)]:
    subs = [frozenset(s) for s in itertools.combinations(range(k), m)]
    total = 0; cnt = {}; greedy_ok = True
    for e in itertools.product(subs, repeat=alpha):
        x = chromatic_index(list(e)); cnt[x] = cnt.get(x, 0) + 1; total += 1
        Delta = max(sum(1 for s in e if v in s) for v in range(k))
        greedy_ok &= x <= m * (Delta - 1) + 1
    tail = {r: F(sum(v for x, v in cnt.items() if x >= r), total) for r in range(1, alpha + 1)}
    lemma = all(tail[r] >= F(m, k) ** r for r in tail)
    print(f"(a) k={k} m={m} alpha={alpha}: distribution of chi' {dict(sorted(cnt.items()))}; "
          f"Pr[chi'>=r] >= (m/k)^r for r=1..{alpha}: {lemma}; greedy bound (77) holds on every sample: {greedy_ok}")
    okall &= lemma and greedy_ok
for L in range(1, 5):
    q = 32 * (L + 1)
    t = 2 ** q
    k = 3 * t * t + 6 * q
    r = L * t * q + 1
    # q/k >= 2^-(2q+3): k <= 8 t^2 = 2^(2q+3) and q >= 1
    bound1 = k <= 2 ** (2 * q + 3)
    bound2 = (2 * q + 3) * r < t * t
    print(f"(b) L={L}, q={q}: k = 3t^2+6q <= 2^(2q+3): {bound1}; (2q+3) r < t^2: {bound2} "
          f"(log2 of the two sides: {((2 * q + 3) * r).bit_length() - 1} vs {2 * q}); r > L t log2 t: {r > L * t * q}")
    okall &= bound1 and bound2 and r > L * t * q
print("ALL_CHECKS", okall)
raise SystemExit(0 if okall else 1)
