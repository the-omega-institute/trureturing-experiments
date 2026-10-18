"""OEIS A387064: a(n) = number of entries in rows 0..n of Pascal's triangle divisible by n.
Conjecture (C. W. Wu): a(p^k) = (p - 1)*p^(k - 1) for p prime, k >= 1.

Checks the published data (n <= 25) and the conjecture for all prime powers p^k <= N,
together with the proof's two facts: no entry of rows m < p^k is divisible by p^k, and
C(p^k, j) is divisible by p^k exactly when p does not divide j.
"""
import sys
from math import comb
from sympy import primerange

N = int(sys.argv[1]) if len(sys.argv) > 1 else 600
DATA = [0, 3, 1, 2, 2, 4, 3, 6, 4, 6, 10, 10, 12, 12, 21, 22, 8, 16, 18, 18, 30, 42, 47, 22, 38, 20]

def a(n):
    if n == 0:
        return 0
    return sum(1 for m in range(n + 1) for j in range(m + 1) if comb(m, j) % n == 0)

ok_data = [a(n) for n in range(len(DATA))] == DATA
ok_conj = ok_facts = True
for p in primerange(2, N + 1):
    q = p
    while q <= N:
        k = 0; t = q
        while t > 1: t //= p; k += 1
        if a(q) != (p - 1) * p ** (k - 1): ok_conj = False
        if any(comb(m, j) % q == 0 for m in range(q) for j in range(m + 1)): ok_facts = False
        if any((comb(q, j) % q == 0) != (j % p != 0) for j in range(1, q)): ok_facts = False
        q *= p
print(f"N={N}: data match: {ok_data}; a(p^k) = (p-1)p^(k-1) for all p^k <= N: {ok_conj}; proof facts: {ok_facts}")
good = ok_data and ok_conj and ok_facts
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
