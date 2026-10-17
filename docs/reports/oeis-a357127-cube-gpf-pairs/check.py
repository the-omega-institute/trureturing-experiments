"""OEIS A357127: a(n) = A081257(n) if A081257(n) > n, otherwise a(n) = 1 (n >= 2), where A081257(n)
is the greatest prime factor of n^3 - 1.

Checked statements:
- Conjecture 1: if a(n) != 1, then a(n) = a(a(n) - n - 1).
- Conjecture 2: if n != m and a(n) = a(m), then a(n) = gcd(n^2+n+1, m^2+m+1) = n + m + 1.
  (Refuted: a(35) = a(61) = 97, while gcd(1261, 3783) = 1261 = 13 * 97.)

For 2 <= n <= N the script computes a(n) by factorization (a(m) for the partner index m = a(n) - n - 1
is computed directly, so the partner may exceed N). It checks Conjecture 1, lists every pair n < m <= N
with a(n) = a(m) != 1 for which Conjecture 2 fails, and checks the proofs' predictions: a(n) = a(m) = p != 1
with n != m forces n + m + 1 = p, and gcd(n^2+n+1, m^2+m+1) = p * gcd((n^2+n+1)/p, m - n).
It also checks the entry's comment that every prime value occurs exactly twice, for primes p <= N.
"""
import sys
from math import gcd
from collections import defaultdict
from sympy import factorint

N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
DATA = [7,13,7,31,43,19,73,13,37,19,157,61,211,241,1,307,1,127,421,463,1,79,601,31,37]
cache = {}
def a(n):
    if n not in cache:
        g = max(factorint(n ** 3 - 1)) if n >= 2 else 0
        cache[n] = g if g > n else 1
    return cache[n]

good = [a(n) for n in range(2, 2 + len(DATA))] == DATA
print("data:", good)
c1 = [n for n in range(2, N + 1) if a(n) != 1 and not (a(n) - n - 1 >= 2 and a(a(n) - n - 1) == a(n))]
print("Conjecture 1 failures n<=%d:" % N, c1[:10], len(c1))
good = good and c1 == []
occ = defaultdict(list)
for n in range(2, N + 1):
    if a(n) != 1: occ[a(n)].append(n)
fails, pred_ok = [], True
for p, ns in occ.items():
    for i in range(len(ns)):
        for j in range(i + 1, len(ns)):
            n, m = ns[i], ns[j]; A, B = n * n + n + 1, m * m + m + 1; g = gcd(A, B)
            if n + m + 1 != p or g != p * gcd(A // p, m - n): pred_ok = False
            if not (p == g == n + m + 1): fails.append((n, m, p, g))
fails.sort()
print("Conjecture 2 failures among pairs n<m<=%d with a(n)=a(m)!=1:" % N, len(fails), "first:", fails[:5])
print("prediction n+m+1 = p and gcd = p*gcd((n^2+n+1)/p, m-n):", pred_ok)
twice = all(len(occ[p]) == 2 for p in occ if p <= N)
print("every prime value p<=%d occurs exactly twice:" % N, twice)
good = good and pred_ok and twice and fails[:1] == [(35, 61, 97, 1261)]
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
