"""OEIS A347285 (De Vlieger): row n starts with n and continues with e_k = floor(log_{p_k}(p_{k-1}^{e_{k-1}}))
while e_k > 0. Conjecture: for n > 0, n - S = 1, where S is the sum of the absolute values of the first
differences of row n.

Usage: python3 check.py N   (rows 1..N, exact integer comparisons of prime powers)
"""
import sys

def primes(m):
    s = [True] * (m + 1); s[0] = s[1] = False
    for i in range(2, int(m ** 0.5) + 1):
        if s[i]:
            for j in range(i * i, m + 1, i): s[j] = False
    return [i for i in range(m + 1) if s[i]]

P = primes(200000)

def row(n):
    r = [n]; k = 0
    while True:
        p, q = P[k], P[k + 1]
        v = p ** r[-1]          # p_{k-1}^{e_{k-1}}
        e = 0
        while q ** (e + 1) <= v: e += 1   # floor(log_q v); q^(e+1) = v is impossible for distinct primes
        if e == 0: return r
        r.append(e); k += 1

DATA = [0,1,2,1,3,1,4,2,1,5,3,2,1,6,3,2,1,7,4,2,1,8,5,3,2,1,9,5,3,2,1,10,6,4]

def main():
    N = int(sys.argv[1])
    flat = [0] + [x for n in range(1, 12) for x in row(n)]
    data = flat[:len(DATA)] == DATA
    bad, last, mono = [], [], True
    for n in range(1, N + 1):
        r = row(n)
        S = sum(abs(r[i + 1] - r[i]) for i in range(len(r) - 1))
        if n - S != 1: bad.append(n)
        if r[-1] != 1: last.append(n)
        if any(r[i + 1] >= r[i] for i in range(len(r) - 1)): mono = False
    print(f"data: {data}; rows 1..{N}: n - S = 1 fails at {bad[:5]}; last term != 1 at {last[:5]}; strictly decreasing: {mono}; longest row length {len(row(N))}")
    ok = data and bad == [] and last == [] and mono
    print("RESULT", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

main()
