"""OEIS A343881 (Kagey): T(n, k) is the least m > k such that k^x * m^y = c^n for some positive integers
c, x, y with x < n and y < n. Conjecture: T(p, k) = A064549(k) = k * Product_{q | k} q for fixed k > 1
and sufficiently large p.

Usage: python3 check.py K P   (corrected formula for 2 <= k <= K at the primes 31 <= p <= P)
"""
import sys
from math import gcd
from functools import reduce

FS = {}
def factor(n):
    if n in FS: return FS[n]
    m, f, d = n, {}, 2
    while d * d <= m:
        while m % d == 0: f[d] = f.get(d, 0) + 1; m //= d
        d += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    FS[n] = f
    return f

def is_prime(n): return n > 1 and factor(n) == {n: 1}

def works(n, k, m):
    """exists 1 <= x, y < n with k^x m^y an n-th power (exponent test, exact)"""
    fk, fm = factor(k), factor(m)
    qs = set(fk) | set(fm)
    for y in range(1, n):
        for x in range(1, n):
            if all((x * fk.get(q, 0) + y * fm.get(q, 0)) % n == 0 for q in qs):
                return (x, y)
    return None

def T(n, k):
    m = k + 1
    while not works(n, k, m): m += 1
    return m

def a064549(k): return k * reduce(lambda u, v: u * v, factor(k).keys(), 1)

def corrected(k):  # k^(1 + 1/g), g = gcd of the exponents of k
    f = factor(k); g = reduce(gcd, f.values()); r = 1
    for q, e in f.items(): r *= q ** (e // g)
    return k * r

TABLE = {2: [4, 8, 12, 9, 20, 24, 28, 18, 16, 40], 3: [8, 4, 9, 16, 25, 36, 49, 27, 24, 80],
         4: [4, 4, 9, 8, 20, 24, 28, 18, 12, 40], 5: [32, 4, 9, 8, 25, 36, 49, 16, 27, 100],
         6: [4, 4, 9, 8, 20, 24, 28, 9, 16, 40], 7: [128, 4, 9, 8, 25, 36, 49, 16, 27, 100],
         8: [4, 4, 9, 8, 20, 24, 28, 16, 12, 40], 9: [8, 4, 9, 8, 25, 36, 49, 16, 24, 80],
         10: [4, 4, 9, 8, 20, 24, 28, 16, 16, 40], 11: [2048, 4, 9, 8, 25, 36, 49, 16, 27, 100]}

def main():
    K, P = int(sys.argv[1]), int(sys.argv[2])
    tab = all([T(n, k) for k in range(1, 11)] == row for n, row in TABLE.items())
    print(f"table rows n = 2..11, k = 1..10 reproduced: {tab}")
    t12 = [(p, T(p, 12)) for p in range(2, 48) if is_prime(p)]
    print(f"T(p, 12) for primes p < 48: {t12}; A064549(12) = {a064549(12)}")
    never = [n for n in range(2, 61) if works(n, 12, 72)]
    print(f"n <= 60 for which m = 72 works at k = 12: {never}")
    primes = [p for p in range(31, P + 1) if is_prime(p)]
    bad = [(p, k) for p in primes for k in range(2, K + 1) if T(p, k) != corrected(k)]
    print(f"T(p, k) = k^(1+1/g(k)) for 2 <= k <= {K}, primes {primes}: {bad == []} {bad[:5]}")
    diff = [k for k in range(2, K + 1) if corrected(k) != a064549(k)]
    print(f"k <= {K} where k^(1+1/g) != A064549(k) (unequal exponents): {diff[:12]}{' ...' if len(diff) > 12 else ''} ({len(diff)})")
    ok = tab and all(t == 144 for p, t in t12 if p >= 7) and never == [] and bad == [] and diff[0] == 12
    print("RESULT", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

main()
