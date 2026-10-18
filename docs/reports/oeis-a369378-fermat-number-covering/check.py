"""OEIS A369378 (T. Ordowski): a(n) is the smallest k such that 2^(2^n+1) - 1 + 2^k is prime.
Conjecture: if n > 0, then 2^(2^n+1) - 1 + 2^k is not prime for every k < 2^n.

Proof checked here (Fermat-number covering): put m = 2^n + 1 and N = 2^m - 1 + 2^k.
- k = 0: N = 2^m is not prime.
- 0 < k < 2^n: let j = v_2(k) < n and F_j = 2^(2^j) + 1. Then 2^(2^n) = 1 (mod F_j), so 2^m = 2,
  and 2^k = (2^(2^j))^(k/2^j) = -1 (mod F_j) because k/2^j is odd. Hence N = 0 (mod F_j), and
  N > F_j, so N is composite.

For 1 <= n <= NMAX the script checks every k < 2^n: that F_{v_2(k)} divides N, that N > F_{v_2(k)}, and
for k = 0 that N = 2^m. It also recomputes the published terms a(1..3), a(5), a(6) (and a(4) = 1968 when
run with --full), and checks that k = 2^n gives a prime for n = 1, 2 (so the bound k < 2^n is sharp there).
"""
import sys
from sympy import isprime
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 12
FULL = '--full' in sys.argv
ok = True
for n in range(1, NMAX + 1):
    m = 2 ** n + 1
    for k in range(0, 2 ** n):
        N = 2 ** m - 1 + 2 ** k
        if k == 0:
            ok = ok and N == 2 ** m
            continue
        j = (k & -k).bit_length() - 1
        F = 2 ** (2 ** j) + 1
        if N % F != 0 or N <= F: ok = False; print('FAIL', n, k)
print(f"covering holds for 1 <= n <= {NMAX}, all k < 2^n: {ok}")
def a(n, kmax):
    m = 2 ** n + 1
    for k in range(kmax + 1):
        if isprime(2 ** m - 1 + 2 ** k): return k
    return None
data = {1: 2, 2: 4, 3: 16, 5: 64, 6: 64}
if FULL: data[4] = 1968
got = {n: a(n, v) for n, v in data.items()}
ok_data = got == data
print("published terms", data, "recomputed", got, ":", ok_data)
sharp = [isprime(2 ** (2 ** n + 1) - 1 + 2 ** (2 ** n)) for n in (1, 2)]
print("k = 2^n prime for n = 1, 2:", sharp)
good = ok and ok_data and all(sharp)
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
