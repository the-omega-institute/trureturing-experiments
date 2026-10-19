#!/usr/bin/env python3
"""Variant of check.py: A(x)^2 = A((x + c A(x)^2)^2) for rational c; checks the defining equation exactly to order 70
and the vanishing [x^(2^n-2^k)] (x/A)^(2^n) = 0 for 2^n <= 64, k <= n-2. Usage: check_c.py c1 c2 ..."""
import sys
from fractions import Fraction as Fr
M = 70
def mul(a, b, m):
    c = [Fr(0)] * (m + 1)
    for i, x in enumerate(a[:m + 1]):
        if x:
            for j in range(0, m + 1 - i):
                if b[j]: c[i + j] += x * b[j]
    return c
def compose(a, g, m):
    res = [Fr(0)] * (m + 1); p = [Fr(1)] + [Fr(0)] * m
    for i in range(1, m + 1):
        p = mul(p, g, m)
        if a[i]:
            for t in range(m + 1): res[t] += a[i] * p[t]
    return res
ok = True
for cs in sys.argv[1:]:
    c = Fr(cs)
    # solve coefficient by coefficient: A = x + sum a_t x^t, unknown a_t enters linearly at order t+1
    A = [Fr(0), Fr(1)] + [Fr(0)] * (M - 1)
    for t in range(2, M + 1):
        def resid(val):
            A[t] = val
            A2 = mul(A, A, t + 1); inner = [(c * A2[s] + (1 if s == 1 else 0)) for s in range(t + 2)]
            return (mul(A, A, t + 1)[t + 1] - compose(A + [Fr(0)], mul(inner, inner, t + 1), t + 1)[t + 1])
        r0 = resid(Fr(0)); r1 = resid(Fr(1)); A[t] = -r0 / (r1 - r0)
    A2 = mul(A, A, M); inner = [(c * A2[s] + (1 if s == 1 else 0)) for s in range(M + 1)]
    eq = mul(A, A, M)[:M] == compose(A, mul(inner, inner, M), M)[:M]
    u = A[1:] + [Fr(0)]; xA = [Fr(0)] * (M + 1); xA[0] = Fr(1)
    for t in range(1, M + 1): xA[t] = -sum(u[i] * xA[t - i] for i in range(1, t + 1))
    bad = []
    for n in range(2, 7):
        N = 2 ** n; Q = [Fr(1)] + [Fr(0)] * N
        for _ in range(N): Q = mul(Q, xA, N)
        bad += [(n, k) for k in range(n - 1) if Q[N - 2 ** k] != 0]
    print(f"c={c}: equation_ok={eq} a2..a4={A[2:5]} failures={bad}")
    ok = ok and eq and not bad
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
