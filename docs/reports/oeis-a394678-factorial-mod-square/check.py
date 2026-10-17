"""OEIS A394678: a(n) = (n! - 1) mod n^2.
Conjecture: excluding a(4) = 7, if a(n) is prime then n is a prime with n^2 - n - 1 prime.

For 1 <= n <= N the script checks the conjecture, the published data, and the two
lemmas of the proof: a(p) = p^2 - p - 1 for prime p (Wilson), and n^2 | n! for
composite n > 4, so a(n) = n^2 - 1.
"""
import sys
from sympy import isprime

N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
DATA = [0, 1, 5, 7, 19, 35, 41, 63, 80, 99, 109, 143, 155, 195, 224, 255, 271, 323, 341, 399]
f, vals, ok_conj, ok_lem = 1, [], True, True
for n in range(1, N + 1):
    f *= n
    a = (f - 1) % (n * n)
    vals.append(a)
    if n != 4 and isprime(a) and not (isprime(n) and isprime(n * n - n - 1)):
        ok_conj = False
    if isprime(n) and a != n * n - n - 1:
        ok_lem = False
    if n > 4 and not isprime(n) and (f % (n * n) or a != n * n - 1):
        ok_lem = False
ok_data = vals[:len(DATA)] == DATA
print(f"N={N}: conjecture holds: {ok_conj}; lemmas hold: {ok_lem}; data match: {ok_data}")
good = ok_conj and ok_lem and ok_data
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
