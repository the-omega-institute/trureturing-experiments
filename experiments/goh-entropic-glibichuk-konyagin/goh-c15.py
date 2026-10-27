"""Exact check of the zero-atom counterexample to Goh (2026) Conjecture 15.

For k and a prime p, mu = (1-1/k) delta_0 + (1/k) U_p.  The script computes the
law of Z = sum_{i<=k} X_i X'_i exactly over the rationals, the Shannon entropies
in bits in double precision, and the margin min(2H(X), log2 p) - 1 - H(Z).
"""
from fractions import Fraction as F
from math import isqrt
from math import log2
import sys


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))

def H(ps):
    return -sum(float(q) * log2(float(q)) for q in ps if q)

def run(k, p):
    assert is_prime(p)
    mu = [F(k - 1, k) + F(1, k * p)] + [F(1, k * p)] * (p - 1)
    prod = [F(0)] * p
    for a in range(p):
        for b in range(p):
            prod[a * b % p] += mu[a] * mu[b]
    z = [F(1)] + [F(0)] * (p - 1)
    for _ in range(k):
        out = [F(0)] * p
        for a in range(p):
            if z[a]:
                for b in range(p):
                    out[(a + b) % p] += z[a] * prod[b]
        z = out
    assert sum(z) == 1
    hx, hz, L = H(mu), H(z), log2(p)
    margin = min(2 * hx, L) - 1 - hz
    return hx, hz, margin

if __name__ == "__main__":
    ok = True
    for k, p in [(3, 67), (3, 257), (4, 257), (5, 1031)]:
        hx, hz, m = run(k, p)
        print(f"k={k} p={p} H(X)={f'{hx:.12f}'} H(Z)={f'{hz:.12f}'} margin={f'{m:.12f}'}")
        ok &= m > 0
    print("ALL_VIOLATE" if ok else "SOME_HOLD")
    sys.exit(0 if ok else 1)
