"""Cullinan-Young, arXiv:2603.05408v1, Remark 8: for even N, the leading coefficient of
I_N(x) = sum_{i=0}^N sgn(i - N/2) l_i(x), l_i(x) = prod_{j != i} (x - (j - N/2)) / (i - j),
is (-1)^((N-2)/2) C((N-2)/2) / (N-1)!, with C(n) the n-th Catalan number ("Numerical experimentation suggests").
The program expands I_N exactly (fractions) for N = 2, 4, ..., 40 and checks that deg I_N = N - 1 and that the
coefficient of x^(N-1) equals the formula. Standard library only.
"""
from fractions import Fraction as Fr
from math import comb, factorial
def polymul(p, q):
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): r[i + j] += a * b
    return r
def I(N):
    nodes = [Fr(j) - Fr(N, 2) for j in range(N + 1)]
    total = [Fr(0)] * (N + 1)
    for i in range(N + 1):
        s = (nodes[i] > 0) - (nodes[i] < 0)
        if s == 0: continue
        p = [Fr(1)]
        for j in range(N + 1):
            if j != i: p = polymul(p, [-nodes[j] / (i - j), Fr(1) / (i - j)])
        total = [a + s * b for a, b in zip(total, p)]
    return total
ok = True
for N in range(2, 41, 2):
    c = I(N); m = N // 2
    deg = max(k for k, v in enumerate(c) if v != 0)
    want = Fr((-1) ** (m - 1) * comb(2 * m - 2, m - 1), m * factorial(N - 1))
    good = deg == N - 1 and c[N - 1] == want
    ok &= good
    if N <= 8 or not good: print("N=%d deg=%d lead=%s formula=%s %s" % (N, deg, c[N - 1], want, "ok" if good else "MISMATCH"))
print("checked N = 2..40 even")
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
