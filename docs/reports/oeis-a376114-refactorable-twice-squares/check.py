"""OEIS A376114 (W. Puszkarz): refactorable numbers (tau(k) | k) that are twice a square.
Conjecture: a(n) = A181795(n)/2, where A181795 lists the k such that the number of odd
divisors of k is an odd divisor of k and the number of even divisors of k is an even
divisor of k.

For k < N the script builds both sequences from their definitions, checks the conjecture
term by term, checks that every term of A181795 is an even square (the proof's first step),
and checks the proof's criterion: for m = 2^a * u (u odd) and c the odd part of 2a + 2,
2m^2 is refactorable iff 4m^2 is in A181795 iff c * tau(u^2) divides u^2.
"""
import sys
from math import isqrt

N = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
tau = [0] * (2 * N + 1)
todd = [0] * (2 * N + 1)
for d in range(1, 2 * N + 1):
    for j in range(d, 2 * N + 1, d):
        tau[j] += 1
        if d % 2: todd[j] += 1

def refactorable(k): return k % tau[k] == 0
def in181795(k):
    o, e = todd[k], tau[k] - todd[k]
    return o % 2 == 1 and k % o == 0 and e > 0 and e % 2 == 0 and k % e == 0

A376114 = [k for k in range(1, N) if refactorable(k) and (k % 2 == 0 and isqrt(k // 2) ** 2 == k // 2)]
A181795 = [k for k in range(1, 2 * N) if in181795(k)]
ok_conj = all(2 * a == b for a, b in zip(A376114, A181795)) and len(A181795) >= len(A376114)
ok_sq = all(k % 2 == 0 and isqrt(k) ** 2 == k for k in A181795)
ok_crit = True
for m in range(1, isqrt(N // 2) + 1):
    a, u = 0, m
    while u % 2 == 0: u //= 2; a += 1
    c = 2 * a + 2
    while c % 2 == 0: c //= 2
    crit = (u * u) % (c * tau[u * u]) == 0
    if refactorable(2 * m * m) != crit or in181795(4 * m * m) != crit: ok_crit = False
print(f"N={N}: {len(A376114)} terms of A376114 below N; a(n) = A181795(n)/2 for all of them: {ok_conj}; "
      f"A181795 terms are even squares: {ok_sq}; criterion c*tau(u^2) | u^2 matches both: {ok_crit}")
print("A376114 head:", A376114[:12])
good = ok_conj and ok_sq and ok_crit and A376114[:16] == [2,8,18,72,128,288,450,882,1152,1250,1800,2178,3042,3528,4050,5000]
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
