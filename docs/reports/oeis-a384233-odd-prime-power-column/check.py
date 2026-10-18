"""OEIS A384233 (O. E. Pol): T(n,k) = n-th number whose largest odd noncomposite divisor
(1 or an odd prime) is its k-th divisor. Conjecture: T(n,2) = A061345(n), i.e. column 2
lists the powers p^e (e >= 1) of odd primes p in increasing order.

The script computes, for every m < M, the position of the largest odd noncomposite
divisor among the sorted divisors of m; checks that position 2 occurs exactly for the odd
prime powers; rebuilds the published antidiagonals of the array; and compares column 2
with A061345 (offset 0: 1, 3, 5, 7, 9, ...).
"""
import sys
from sympy import divisors, isprime, factorint

M = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
DATA = [1,2,3,4,5,6,8,7,10,20,16,9,12,28,42,32,11,14,30,60,84,64,13,15,40,66]
A061345 = [1,3,5,7,9,11,13,17,19,23,25,27,29,31,37,41,43,47,49,53,59,61]

def pos(m):
    ds = divisors(m)
    L = max(d for d in ds if d % 2 == 1 and (d == 1 or isprime(d)))
    return ds.index(L) + 1

def odd_prime_power(m):
    f = factorint(m)
    return m > 1 and len(f) == 1 and 2 not in f

cols = {}
ok_char = True
for m in range(1, M):
    k = pos(m)
    cols.setdefault(k, []).append(m)
    if (k == 2) != odd_prime_power(m): ok_char = False
T = lambda n, k: cols[k][n - 1]
anti = []
for s in range(2, 9):            # upward antidiagonals n + k = s: T(s-1,1), ..., T(1,s-1)
    for k in range(1, s):
        anti.append(T(s - k, k))
ok_data = anti[:len(DATA)] == DATA
ok_col = cols[2][:len(A061345) - 1] == A061345[1:]
print(f"M={M}: position 2 <=> odd prime power for all m < M: {ok_char}; antidiagonals match data: {ok_data}; column 2 = A061345(n), n >= 1: {ok_col}")
good = ok_char and ok_data and ok_col
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
