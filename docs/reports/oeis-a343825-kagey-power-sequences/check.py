"""OEIS A343825 (Kagey): T(n, k) is the least m such that there is a sequence k = b_1 <= ... <= b_t = m in
which no term appears n or more times and whose product is an n-th power. Conjecture: T(p, k) = A071537(k)
for fixed k and sufficiently large prime p.

For a prime p and k >= 2, a sequence ending at m is a multiplicity vector c_j in [0, p-1] (j = k..m) with
c_k, c_m >= 1 and sum c_j v(j) = 0 mod p, v(j) the exponent vector of j. Such a vector exists iff v(k) lies in
the F_p-span of v(k+1..m) and v(m) in the F_p-span of v(k..m-1) (a vector space is not the union of two
proper subspaces). A071537(k) is the least m with v(k) in the Q-span of v(k+1..m).

Usage: python3 check.py K   (compare for 1 <= k <= K at the primes 31, 37, 41, 43)
"""
import sys
from fractions import Fraction

FS = {}
def factor(n):
    if n in FS: return FS[n]
    m, f, d = n, {}, 2
    while d * d <= m:
        while m % d == 0: f[d] = f.get(d, 0) + 1; m //= d
        d += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    FS[n] = f
    return f

def in_span(vecs, v, p=None):
    """v in the span of vecs over Q (p is None) or F_p; vectors are dicts prime -> exponent"""
    def norm(u):
        return {q: (e if p is None else e % p) for q, e in u.items() if (e if p is None else e % p)}
    basis = []
    def red(u):
        u = dict(u)
        for piv, b in basis:
            if u.get(piv, 0):
                f = Fraction(u[piv]) / b[piv] if p is None else u[piv] * pow(b[piv], -1, p) % p
                for q, x in b.items():
                    nv = u.get(q, 0) - f * x
                    if p is not None: nv %= p
                    if nv: u[q] = nv
                    else: u.pop(q, None)
        return u
    for w in vecs:
        u = red(norm(w))
        if u: basis.append((min(u), u))
    return not red(norm(v))

def A071537(k):
    if k == 1: return 1
    m = k
    while True:
        m += 1
        if in_span([factor(j) for j in range(k + 1, m + 1)], factor(k)): return m

def T(p, k):
    if k <= 1: return k
    if all(e % p == 0 for e in factor(k).values()): return k
    m = k
    while True:
        m += 1
        V = [factor(j) for j in range(k, m + 1)]
        if in_span(V[1:], V[0], p) and in_span(V[:-1], V[-1], p): return m

A071537_DATA = [1,4,6,8,10,9,14,12,15,16,22,18,26,21,20,24,34,25,38,27,28,33,46,30,32,39,35,36,58,40,62,42,44,51,45,48,74,57,52,49,82,50,86,55,54,69,94,56,60,63,68,65,106,64,66,70,76,87,118,72,122,93,75,77,78,80,134]
ROWS = {2: [0,1,6,8,4,10,12,14,15,9,18], 3: [0,1,4,6,9,10,12,14,8,16,15], 5: [0,1,4,6,8,10,9,14,12,15,16], 7: [0,1,4,6,8,10,9,14,12,15,16]}

def main():
    K = int(sys.argv[1])
    d1 = [A071537(k) for k in range(1, len(A071537_DATA) + 1)] == A071537_DATA
    d2 = all([T(p, k) for k in range(0, 11)] == row for p, row in ROWS.items())
    print(f"A071537 data: {d1}; A343825 prime rows 2, 3, 5, 7 (k = 0..10): {d2}")
    res = {}
    for p in [31, 37, 41, 43]:
        res[p] = [k for k in range(1, K + 1) if T(p, k) != A071537(k)]
        print(f"p = {p}: k <= {K} with T(p, k) != A071537(k): {res[p]}")
    small = {k: [p for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31] if T(p, k) != A071537(k)] for k in [2, 4, 8, 9, 12, 16, 27, 32, 64]}
    print(f"primes p <= 31 with T(p, k) != A071537(k): {small}")
    # the exceptions at p = 37 disappear at the larger primes 41 and 43
    ok = d1 and d2 and res[31] == [] and res[41] == [] and res[43] == [] and all(T(41, k) == A071537(k) == T(43, k) for k in res[37])
    print("RESULT", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

main()
