"""Exact check of the rank-sequence family s_p(n) = p*8^n - 4^n - 2^n (arXiv:2501.14643v1, Section 1 conjecture).

s_p(n)^M = (-1)^M sum_j a_{p,M,j} (2^{M+j})^n with (1 + x - p x^2)^M = sum_j a_{p,M,j} x^j, so the rank of s_p^M
is the number of nonzero a_{p,M,j}.  The script checks this against the minimal recurrence order computed by
Berlekamp-Massey over Q on the integer terms, and checks that the rank sequences of the primes
3, 11, 29, 61, 127 are pairwise distinct.
"""
from fractions import Fraction as F


def poly_pow(c, M):
    r = [1]
    for _ in range(M):
        n = [0] * (len(r) + len(c) - 1)
        for i, a in enumerate(r):
            for j, b in enumerate(c):
                n[i + j] += a * b
        r = n
    return r


def support_rank(p, M):
    return sum(1 for a in poly_pow([1, 1, -p], M) if a != 0)


def bm_order(seq):
    s = [F(x) for x in seq]
    C, B, L, m, b = [F(1)], [F(1)], 0, 1, F(1)
    for n in range(len(s)):
        d = s[n] + sum(C[i] * s[n - i] for i in range(1, L + 1))
        if d == 0:
            m += 1
            continue
        T = C[:]
        C = C + [F(0)] * (len(B) + m - len(C))
        for i in range(len(B)):
            C[i + m] -= d / b * B[i]
        if 2 * L <= n:
            L, B, b, m = n + 1 - L, T, d, 1
        else:
            m += 1
    return L


def s(p, n):
    return p * 8 ** n - 4 ** n - 2 ** n


ok = True
for p in [3, 5, 7, 11, 29]:
    terms = [s(p, n) for n in range(8)]
    rec = all(s(p, n + 3) == 14 * s(p, n + 2) - 56 * s(p, n + 1) + 64 * s(p, n) for n in range(40))
    pos = all(s(p, n) > 0 for n in range(40))
    print(f"p={p}: s_p starts {terms[:4]}, recurrence(14,-56,64) {rec}, positive {pos}, rank(s_p)={bm_order([s(p, n) for n in range(20)])}")
    ok &= rec and pos
for p in [3, 5, 7, 11]:
    seq = [support_rank(p, M) for M in range(1, 2 * p + 3)]
    first_drop = next(M for M in range(1, 2 * p + 3) if support_rank(p, M) < 2 * M + 1)
    print(f"p={p}: R_p(M) for M=1..{2 * p + 2}: {seq}; first M with R_p(M) < 2M+1: {first_drop}")
    ok &= all(support_rank(p, M) == 2 * M + 1 for M in range(1, p))
    ok &= support_rank(p, 2 * p + 1) <= 2 * (2 * p + 1) - 1
for p, M in [(3, 6), (3, 7), (5, 4), (5, 11), (7, 15), (11, 7)]:
    terms = [s(p, n) ** M for n in range(4 * (2 * M + 1) + 4)]
    a, b = support_rank(p, M), bm_order(terms)
    print(f"p={p} M={M}: coefficient support {a}, Berlekamp-Massey order {b}, 2M+1={2 * M + 1}")
    ok &= a == b
print("coefficients of (1+x-3x^2)^7:", poly_pow([1, 1, -3], 7))


def isprime(n):
    return n > 1 and all(n % k for k in range(2, int(n ** 0.5) + 1))


chain = [3]
while len(chain) < 5:
    q = 2 * chain[-1] + 2
    while not isprime(q):
        q += 1
    chain.append(chain[-1] * 0 + q)
print("prime chain p_{i+1} = least prime > 2 p_i + 1:", chain)
for i, pi in enumerate(chain):
    for pj in chain[i + 1:]:
        M = 2 * pi + 1
        ri, rj = support_rank(pi, M), support_rank(pj, M)
        print(f"  M={M}: R_{pi}(M)={ri}, R_{pj}(M)={rj}")
        ok &= ri < rj
print("ALL_CHECKS", ok)
raise SystemExit(0 if ok else 1)
