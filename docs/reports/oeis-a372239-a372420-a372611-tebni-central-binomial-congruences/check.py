"""OEIS A372239, A372420, A372611: a_c(n) = [x^n] (1 + c x)/((1 - 2x) sqrt(1 - 4x)) for c = 2, 1, 3.
Conjectures: for primes p == 1 (mod 4), a_c(p) == 4 + c (mod p); for p == 3 (mod 4),
a_c(p) == -c (mod p).  (c = 2: 6 and -2; c = 1: 5 and -1; c = 3: 7 and -3.)

The script computes a_c(n) = b(n) + c*b(n-1), b(n) = Sum_{k<=n} C(2k,k) 2^(n-k), checks the published
data, the conjectures for all primes p < P, and the proof's two congruences:
b(p) == 2*((-1)^((p-1)/2) + 1) and b(p-1) == (-1)^((p-1)/2) (mod p).
"""
import sys
from math import comb
from sympy import primerange
P = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
DATA = {2: [1,6,22,76,262,916,3260,11800,43334,161028,604052,2283048,8681116],
        1: [1,5,18,62,214,750,2676,9708,35718,132926,499228,1888644,7186876],
        3: [1,7,26,90,310,1082,3844,13892,50950,189130,708876,2677452,10175356]}
def b_exact(n): return sum(comb(2*k, k) * 2**(n - k) for k in range(n + 1))
ok = True
for c, data in DATA.items():
    vals = [b_exact(n) + c * (b_exact(n - 1) if n else 0) for n in range(len(data))]
    if vals != data: ok = False; print("data", c)
okc = okl = True
for p in primerange(3, P):
    # b(p) and b(p-1) modulo p
    bp = sum(comb(2*k, k) % p * pow(2, p - k, p) for k in range(p + 1)) % p
    bq = sum(comb(2*k, k) % p * pow(2, p - 1 - k, p) for k in range(p)) % p
    eps = 1 if p % 4 == 1 else -1
    if bp != (2 * (eps + 1)) % p or bq != eps % p: okl = False
    for c in (1, 2, 3):
        want = (4 + c) % p if p % 4 == 1 else (-c) % p
        if (bp + c * bq) % p != want: okc = False
print(f"data match: {ok}; conjectures for 3 <= p < {P}: {okc}; proof congruences: {okl}")
good = ok and okc and okl
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
