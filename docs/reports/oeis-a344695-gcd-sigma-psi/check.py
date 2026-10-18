"""OEIS A344695: a(n) = gcd(sigma(n), psi(n)), psi the Dedekind psi function.
Karttunen-Munn conjecture: the set of primes that appear in the sequence is the set of odd primes.

Usage: python3 check.py N Q
  N: sieve bound for a(n) != 2 and the data
  Q: odd primes q <= Q for the construction n = q^2 r^(q-1)
"""
import sys
from math import gcd

DATA = [1,3,4,1,6,12,8,3,1,18,12,4,14,24,24,1,18,3,20,6,32,36,24,12,1,42,4,8,30,72,32,3,48,54,48,1,38,60,56,18,42,96,44,12,6,72,48,4,1,3,72,14,54,12,72,24,80,90,60,24,62,96,8,1,84,144,68,18,96,144,72,3,74,114,4,20,96,168,80,6,1,126,84,32,108]

def is_prime(m):
    if m < 2: return False
    if m % 2 == 0: return m == 2
    d = 3
    while d * d <= m:
        if m % d == 0: return False
        d += 2
    return True

def sigma_psi_of_factorization(f):
    s = p = 1
    for q, e in f.items():
        s *= (q ** (e + 1) - 1) // (q - 1)
        p *= q ** (e - 1) * (q + 1)
    return s, p

def main():
    N, Q = int(sys.argv[1]), int(sys.argv[2])
    spf = list(range(N + 1))
    for i in range(2, int(N ** 0.5) + 1):
        if spf[i] == i:
            for m in range(i * i, N + 1, i):
                if spf[m] == m: spf[m] = i
    def fac(n):
        f = {}
        while n > 1:
            q = spf[n]; f[q] = f.get(q, 0) + 1; n //= q
        return f
    a = [0] + [gcd(*sigma_psi_of_factorization(fac(n))) for n in range(1, N + 1)]
    data = a[1:len(DATA) + 1] == DATA
    two = [n for n in range(1, N + 1) if a[n] == 2]
    first = {}
    for n in range(1, N + 1):
        if is_prime(a[n]) and a[n] not in first: first[a[n]] = n
    print(f"data: {data}; n <= {N} with a(n) = 2: {two[:5]}; first n with a(n) = q for primes q <= 23: { {q: first.get(q) for q in [2,3,5,7,11,13,17,19,23]} }")
    # the dichotomy behind "2 never occurs": square n has odd sigma; otherwise a prime p with odd exponent gives (p+1) | a(n)
    dich = all((a[n] % 2 == 1) if all(e % 2 == 0 for e in fac(n).values()) else any(a[n] % (p + 1) == 0 for p, e in fac(n).items() if e % 2) for n in range(1, N + 1))
    print(f"square n: a(n) odd; non-square n: (p+1) | a(n) for a prime p of odd exponent (n <= {N}): {dich}")
    # construction: q odd prime, r the least prime = 1 mod q(q+1)(1+q+q^2) with r > 1+q+q^2, n = q^2 r^(q-1)
    cons, bad = [], []
    for q in range(3, Q + 1):
        if not is_prime(q): continue
        M = q * (q + 1) * (1 + q + q * q); r = M + 1
        while not (is_prime(r) and r > 1 + q + q * q): r += M
        s, p = sigma_psi_of_factorization({q: 2, r: q - 1})
        g = gcd(s, p)
        cons.append((q, r)); 
        if g != q: bad.append((q, r, g))
    print(f"construction gcd(sigma, psi)(q^2 r^(q-1)) = q for the {len(cons)} odd primes q <= {Q}: {bad == []}; first pairs (q, r): {cons[:6]}")
    ok = data and two == [] and dich and bad == []
    print("RESULT", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

main()
