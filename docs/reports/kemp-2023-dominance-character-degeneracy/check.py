#!/usr/bin/env python3
"""Exact checks for Kemp, arXiv:2305.06768v2, Conjecture 1 (eq. (3.7)): k_deg <= k_swap + 1.

1. Character bridge: for n <= 10, normalized characters chi_Y(T_k)/dim Y (Murnaghan-Nakayama, exact)
   are degenerate for k = 2..m iff the content power sums <c^q> agree for q = 1..m-1 (eq. (5.22)).
2. Conjecture: for every pair of distinct partitions of n <= N_MAX, the first q with <c^q>(X) != <c^q>(Y)
   satisfies q <= k_swap + 1, i.e. k_deg <= k_swap + 1.
3. Stronger statement: for random decreasing integer vectors a != b with equal length and sum, some power
   sum of degree j in {2, ..., s + 2} differs, s = sign changes of the prefix-sum differences.
Prints ALL_OK on success; exit code 0.
"""
import sys, random
from fractions import Fraction
from functools import lru_cache
from math import factorial

N_MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k): yield (k,) + rest

def contents(lam):
    return [j - i for i, r in enumerate(lam) for j in range(r)]

def content_sums(lam, qmax):
    c = contents(lam)
    return [sum(x ** q for x in c) for q in range(qmax + 1)]

def kswap(x, y):
    d = max(len(x), len(y)); x = list(x) + [0] * (d - len(x)); y = list(y) + [0] * (d - len(y))
    sx = sy = 0; signs = []
    for i in range(d):
        sx += x[i]; sy += y[i]
        if sx != sy: signs.append(1 if sx > sy else -1)
    return sum(1 for a, b in zip(signs, signs[1:]) if a != b)

# --- Murnaghan-Nakayama (exact) via border strips on beta-sets ---
@lru_cache(maxsize=None)
def mn_char(lam, mu):
    """chi^lam(mu): lam, mu partitions (tuples) of the same n."""
    if not mu: return 1 if not lam else 0
    k, rest = mu[0], mu[1:]
    L = len(lam); beta = [lam[i] + (L - 1 - i) for i in range(L)]
    bs = set(beta); total = 0
    for b in beta:
        if b - k >= 0 and (b - k) not in bs:
            height = sum(1 for c in beta if b - k < c < b)
            nb = sorted([c for c in beta if c != b] + [b - k], reverse=True)
            m = len(nb); new = tuple(v for v in (nb[i] - (m - 1 - i) for i in range(m)) if v > 0)
            total += (-1) ** height * mn_char(new, rest)
    return total

def class_size(n, k):  # number of k-cycles (cycle type (k,1^{n-k}))
    return factorial(n) // (factorial(n - k) * k)

def normalized_char(lam, k):
    n = sum(lam)
    mu = tuple([k] + [1] * (n - k)) if k > 1 else tuple([1] * n)
    return Fraction(class_size(n, k) * mn_char(lam, mu), mn_char(lam, tuple([1] * n)))

def check_bridge(nmax=10):
    for n in range(2, nmax + 1):
        P = list(partitions(n))
        chars = {p: [normalized_char(p, k) for k in range(2, n + 1)] for p in P}
        cs = {p: content_sums(p, n) for p in P}
        for i, x in enumerate(P):
            for y in P[i + 1:]:
                for m in range(2, n + 1):
                    a = all(chars[x][k - 2] == chars[y][k - 2] for k in range(2, m + 1))
                    b = all(cs[x][q] == cs[y][q] for q in range(1, m))
                    if a != b:
                        print('BRIDGE_FAIL', n, x, y, m); return False
    print('bridge ok: n <= %d, chars T_2..T_m degenerate iff <c^1..c^(m-1)> agree' % nmax)
    return True

def check_conjecture(nmax):
    worst = {}
    for n in range(2, nmax + 1):
        P = list(partitions(n)); qmax = n + 1
        cs = {p: content_sums(p, qmax) for p in P}
        for i, x in enumerate(P):
            for y in P[i + 1:]:
                q = next((q for q in range(1, qmax + 1) if cs[x][q] != cs[y][q]), None)
                if q is None: print('NO_SEPARATION', x, y); return False
                kdeg = q   # T_2..T_q degenerate, T_{q+1} distinct
                s = kswap(x, y)
                if kdeg > s + 1: print('CONJECTURE_FAIL', n, x, y, 'kdeg', kdeg, 'kswap', s); return False
                worst[s] = max(worst.get(s, 0), kdeg - (s + 1))
        print('n=%d pairs=%d ok' % (n, len(P) * (len(P) - 1) // 2), flush=True)
    print('conjecture ok for all pairs, n <= %d; max(kdeg-(kswap+1)) by kswap:' % nmax, dict(sorted(worst.items())))
    return True

def check_vectors(trials=20000, seed=1):
    rng = random.Random(seed)
    for _ in range(trials):
        d = rng.randint(2, 9)
        a = sorted([rng.randint(-12, 12) for _ in range(d)], reverse=True)
        b = sorted([rng.randint(-12, 12) for _ in range(d)], reverse=True)
        diff = sum(a) - sum(b)
        if diff: b[-1] += diff; b.sort(reverse=True)
        if a == b: continue
        A = B = 0; signs = []
        for i in range(d - 1):
            A += a[i]; B += b[i]
            if A != B: signs.append(1 if A > B else -1)
        s = sum(1 for u, v in zip(signs, signs[1:]) if u != v)
        if not any(sum(x ** j for x in a) != sum(x ** j for x in b) for j in range(2, s + 3)):
            print('VECTOR_FAIL', a, b, s); return False
    print('stronger vector statement ok on %d random pairs' % trials)
    return True

if __name__ == '__main__':
    ok = check_bridge(10) and check_vectors() and check_conjecture(N_MAX)
    print('ALL_OK' if ok else 'FAIL'); sys.exit(0 if ok else 1)
