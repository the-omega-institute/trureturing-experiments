#!/usr/bin/env python3
"""Exact check of the proposed proof of de Oliveira-Barbosa-Galvao, Quantum 8, 1312 (2024),
arXiv:2212.03668v4, Sec. 4.2, Conjecture 2 (symmetric NMQC-xor GHZ resource for the complete
symmetric function C^k is Omega(n^(floor(k/2)-1))).

By Definitions 15-16 and Eq. (72), a symmetric measurement assignment with layer set T computes
C^k_n iff some phase lift g(x) = a_0 + sum_{i in T} a_i * sum_{|S|=i} parity_S(x) is integer-valued
with g(x) = C^k(x) mod 2; it uses q = sum_{i in T} binom(n, i) qubits. Every such g depends on
x only through the weight w, and parity_S = (1 - chi_S)/2 with E_j(w) = sum_{|S|=j} chi_S(x).

Checks, for all 2 <= k <= n <= 24 (m = k // 2, r = k - 2m):
 1. Krawtchouk form: E_j(w) = sum_t (-2)^t binom(n-t, j-t) binom(w, t) equals the direct
    sum_t (-1)^t binom(w, t) binom(n-w, j-t), for all j, w.
 2. C^k(x) = binom(w, k) mod 2 (the parity of k-subsets of the support).
 3. The certificate d (d_{2u+r} = (-1)^(m-u) binom(m, u), u = 0..m, else 0) annihilates the constant
    and every E_j with j < m or j > n - m, and pairs to an odd integer with every integer lift of
    C^k (all lifts agree mod 2 on the support of d, where C^k is 0 for u < m and 1 for u = m).
 4. Hence every symmetric assignment uses a layer j in [m, n - m], so q >= binom(n, m).
"""
import sys
from math import comb

def E_direct(n, j, w):
    return sum((-1)**t * comb(w, t) * comb(n - w, j - t) for t in range(0, j + 1))

def E_kraw(n, j, w):
    return sum((-2)**t * comb(n - t, j - t) * comb(w, t) for t in range(0, j + 1))

ok = True
pairs = 0
for n in range(2, 25):
    for j in range(0, n + 1):
        for w in range(0, n + 1):
            if E_direct(n, j, w) != E_kraw(n, j, w):
                ok = False
    for k in range(2, n + 1):
        m, r = k // 2, k - 2 * (k // 2)
        d = [0] * (n + 1)
        for u in range(0, m + 1):
            d[2*u + r] = (-1)**(m - u) * comb(m, u)
        # annihilation of the constant and of low/high layers
        if sum(d) != 0: ok = False
        for j in range(0, n + 1):
            if j < m or j > n - m:
                if sum(d[w] * E_direct(n, j, w) for w in range(n + 1)) != 0:
                    ok = False
        # pairing with the parity pattern: C^k(w) = binom(w,k) mod 2 on the support
        pattern = [comb(w, k) % 2 for w in range(n + 1)]
        for u in range(0, m):
            if pattern[2*u + r] != 0: ok = False
        if pattern[k] != 1: ok = False
        # sum d_w * g(w) = sum over support of (+-binom) * (even or odd) is odd since binom(m,m) = 1
        pair_parity = sum(abs(d[w]) * pattern[w] for w in range(n + 1)) % 2
        if pair_parity != 1: ok = False
        pairs += 1
print(f"checked {pairs} pairs (k, n) with 2 <= k <= n <= 24: Krawtchouk form, annihilation of the "
      f"constant and of layers j < m or j > n - m, and odd pairing with every lift of C^k")
print("ALL_OK" if ok else "SOME_FAIL")
sys.exit(0 if ok else 1)
