"""Check of the explicit MSTD family for Mandal-Neetu (arXiv:2602.09073) Conjecture 6.2.

Dic_{4n} = <a,b | a^{2n}=1, b^2=a^n, b^{-1}ab=a^{-1}>, elements (i,e) <-> a^i b^e.
For m >= 2, r = m-1 and n > 2r, A_t = {a,...,a^r} u {a^t b}, t in Z/2n.
Claims: |A_t A_t| = 4r, |A_t A_t^{-1}| = 4r-1, so S_n(m) >= 2n.
Also brute-forces S_n(m) for small (n, m) as an independent check of the bound.
"""
from itertools import combinations
import sys

def mul(x, y, n):
    (i, e), (j, f) = x, y
    return ((i + (j if e == 0 else -j) + n * e * f) % (2 * n), (e + f) % 2)

def inv(x, n):
    i, e = x
    return ((-i) % (2 * n), 0) if e == 0 else ((i + n) % (2 * n), 1)

def sizes(A, n):
    P = {mul(x, y, n) for x in A for y in A}
    Q = {mul(x, inv(y, n), n) for x in A for y in A}
    return len(P), len(Q)

ok = True
for m in range(2, 8):
    r = m - 1
    for n in range(2 * r + 1, 2 * r + 9):
        for t in range(2 * n):
            A = [(i, 0) for i in range(1, r + 1)] + [(t, 1)]
            if sizes(A, n) != (4 * r, 4 * r - 1):
                ok = False
                print("FAIL", m, n, t, sizes(A, n))
print("family checked m=2..7, n=2m-1..2m+6:", ok)

G = lambda n: [(i, e) for e in (0, 1) for i in range(2 * n)]
for m, n in [(2, 3), (2, 4), (3, 5), (3, 6), (4, 7)]:
    s = sum(1 for A in combinations(G(n), m) if (lambda p: p[0] > p[1])(sizes(A, n)))
    print(f"S_n(m) brute force m={m} n={n}: {s} >= 2n={2*n}: {s >= 2*n}")
    ok &= s >= 2 * n
print("ALL_OK" if ok else "FAILED")
sys.exit(0 if ok else 1)
