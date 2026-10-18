"""OEIS A073759: "Largest number < n that is neither a divisor of nor relatively prime to n, or 0 if no
such number exists." Conjecture (R. J. Mathar): "a(n) = 2*A046667(n) if n>=5", where A046667(n) = A046666(n)/2
and A046666(n) = n - (smallest prime dividing n).

Route: for composite n >= 5 with least prime factor p, m = n - p is unrelated to n (gcd(m, n) = p < m), and no
k with n - p < k < n shares a prime with n, since such a prime would divide n - k < p. For prime n the set is
empty. n - p is even, so n - p = 2 * ((n - p) / 2).
Checks for 5 <= n <= N: the brute-force a(n) equals 2 * A046667(n); the first OEIS terms agree; a(4) = 0 while
4 - 2 = 2, so the bound n >= 5 is needed.
"""
import math
N = 20000
OEIS = [0,0,0,0,0,4,0,6,6,8,0,10,0,12,12,14,0,16,0,18,18,20,0,22,20,24,24,26,0,28,0,30,30,32,30,34,0,36,36,38,0,40,0,42,42,44,0,46,42,48]

def lpf(n):
    p = 2
    while p * p <= n:
        if n % p == 0:
            return p
        p += 1
    return n

def a(n):
    for k in range(n - 1, 0, -1):
        if n % k and math.gcd(n, k) != 1:
            return k
    return 0

ok = [a(n) for n in range(1, len(OEIS) + 1)] == OEIS
bad = [n for n in range(5, N + 1) if a(n) != 2 * ((n - lpf(n)) // 2)]
print("OEIS terms match:", ok)
print("checked 5..%d, mismatches: %d" % (N, len(bad)))
print("a(4) = %d, 4 - lpf(4) = %d" % (a(4), 4 - lpf(4)))
ok = ok and not bad and a(4) == 0
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
