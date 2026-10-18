"""Hanna OEIS parity/integrality conjectures settled by the research answer (GPT Pro, checked by the orchestrator).

Every sequence is recomputed from its literal definition (or the stated exact rewrite of a bilateral sum),
compared with the OEIS data saved beside this file, and tested against the conjecture:
A393173, A393850 (n >= 1), A382122 (both sign branches), A378258, A375449, A375452, A365574 (Conjecture 1),
A380711, A379764 (both conjectures), and A380681-A380686 (integer iff p prime, with the first non-integer
coefficient for composite p), plus the identities Psi_1 = 0 and Psi_p = 0 (mod p) for random Y.
"""
from fractions import Fraction as Fr
from math import comb, factorial, prod, isqrt
import random

def mul(a, b, n):
    c = [Fr(0)] * (n + 1)
    for i, x in enumerate(a[:n + 1]):
        if x:
            for j, y in enumerate(b[:n + 1 - i]):
                if y: c[i + j] += x * y
    return c
def power(a, k, n):
    r = [Fr(1)] + [Fr(0)] * n
    for _ in range(k): r = mul(r, a, n)
    return r
def inv(a, n):
    b = [Fr(0)] * (n + 1); b[0] = 1 / a[0]
    for k in range(1, n + 1): b[k] = -sum(a[j] * b[k - j] for j in range(1, k + 1)) / a[0]
    return b
def compose(f, g, n):
    res = [Fr(f[0])] + [Fr(0)] * n; P = [Fr(1)] + [Fr(0)] * n
    for k in range(1, n + 1):
        P = mul(P, g, n)
        if k < len(f) and f[k]:
            for i in range(n + 1): res[i] += f[k] * P[i]
    return res
def oeis(name, k):
    v = []
    for line in open(f'oeis-{name}.txt'):
        if line[:2] in ('%S', '%T', '%U'): v += [int(t) for t in line.split(None, 2)[2].strip().strip(',').split(',') if t]
    return v[:k]
def linsolve(A, n, val):
    A[n] = Fr(0); v0 = val(A); A[n] = Fr(1); v1 = val(A); A[n] = -v0 / (v1 - v0)
pow2 = lambda n: n > 0 and n & (n - 1) == 0
integral = lambda A: all(a.denominator == 1 for a in A)
ok = True
def report(name, A, offset, cond, rng):
    global ok
    ref = oeis(name, 12); got = [int(a) for a in A[offset:offset + len(ref)]]
    m = got == ref; c = all(cond(n) for n in rng); i = integral(A)
    print(f'{name}: data match {m}, integral {i}, conjecture on {rng.start}..{rng.stop - 1}: {c}'); ok &= m and c and i

N = 40
A = [Fr(1)] + [Fr(0)] * N                                   # A393173
for n in range(1, N + 1):
    def val(B, n=n):
        d = B[:n + 1]; d[1] -= n * n; return inv(d, n)[n]
    linsolve(A, n, val)
report('a393173', A, 0, lambda n: (int(A[n]) % 2 == 1) == pow2(n), range(1, N + 1))

A = [Fr(1), Fr(1)] + [Fr(0)] * (N - 1)                      # A393850
for n in range(2, N + 1):
    def val(B, n=n):
        P = power(B[:n + 1], n + 1, n); return P[n] - n * (n + 1) * P[n - 1]
    linsolve(A, n, val)
report('a393850', A, 0, lambda n: (int(A[n]) % 2 == 1) == pow2(n), range(1, N + 1))
print('  a393850: a(0) = 1 is odd and 0 is not a power of 2 (the literal reading fails only at n = 0)')

Cat = [Fr(comb(2 * k, k), k + 1) for k in range(N + 2)]      # A382122, both sign branches
for branch, sg in (('oeis', -1), ('other', 1)):
    U = [Fr(1)] + [Fr(0)] * N
    for n in range(1, N + 1):
        s = sum(abs(power(U[:n], m, n)[n + 1 - m]) for m in range(2, n + 2))
        U[n] = sg * (Cat[n + 1] - s)
    A = inv(U, N)
    if branch == 'oeis': report('a382122', A, 0, lambda n: (int(A[n]) % 2 == 1) == pow2(n), range(1, N + 1))
    else:
        c = integral(A) and all((int(A[n]) % 2 == 1) == pow2(n) for n in range(1, N + 1))
        print(f'a382122 other solution branch: integral and same parity {c}'); ok &= c

A = [Fr(0), Fr(1)] + [Fr(0)] * (N - 1)                      # A378258
for n in range(2, N + 1):
    def val(C, n=n):
        C2 = mul(C, C, N); Bs = [C[i] - 5 * C2[i] for i in range(N + 1)]; inner = [C[i] - C2[i] for i in range(N + 1)]
        return compose(Bs, inner, n)[n]
    linsolve(A, n, val)
report('a378258', A, 1, lambda n: (int(A[n]) % 2 == 1) == pow2(n), range(1, N + 1))

def diag(tf, N):                                            # A375449, A375452
    A = [Fr(0)] * (N + 1)
    for n in range(1, N + 1):
        def val(B, n=n):
            xd = [Fr(i) * B[i] for i in range(n + 1)]; den = [Fr(1)] + [-n * B[i] for i in range(1, n + 1)]
            return mul(xd, inv(den, n), n)[n] - n * tf(n)
        linsolve(A, n, val)
    return A
spec = lambda n: any(1 < k and (n == 2 ** k - 2 or n == 2 ** k - 1) for k in range(2, 12))
for name, tf in (('a375449', lambda n: Fr(n) ** n), ('a375452', lambda n: Fr(n) ** (n - 1))):
    A = diag(tf, 60); report(name, A, 1, lambda n: (int(A[n]) % 2 == 1) == spec(n), range(2, 61))

r = [Fr(0)] + [Fr((2 - k) + (k + 1) * (-2) ** (k - 1), 3) for k in range(1, 62)]   # A365574
H = [Fr(0)] * 62; H[1] = Fr(1)
for m in range(2, 62): H[m] = -compose(r, H, m)[m]
A = H[1:]
defok = True
for n in range(2, 20):
    t = [Fr(1)] + [Fr(0)] * n
    for i, v in enumerate([Fr(0)] + A[:n]): t[i] += (n + 1) * v
    defok &= mul(power(t, n, n), power(inv(A[:n + 1], n), n, n), n)[n - 1] == n * Fr(n + 2) ** (n - 2)
print('a365574: explicit reversion satisfies the definition for n = 2..19:', defok); ok &= defok and integral(r)
report('a365574', A, 0, lambda k: (int(A[k]) % 2 == 1) == any(k == 2 ** m - 2 for m in range(1, 12)), range(0, 60))

A = [Fr(1)] + [Fr(0)] * N                                   # A380711
for n in range(1, N + 1):
    Ai = inv(A[:n], n - 1); A3 = mul(mul(Ai, Ai, n - 1), Ai, n - 1)
    A[n] = mul(A[:n], [abs(c) for c in A3], n - 1)[n - 1]
quot = lambda n: Fr(comb(3 * n - 1, n), 3 * n - 1) if n > 0 else Fr(-1)
report('a380711', A, 0, lambda n: quot(n).denominator == 1 and (int(A[n]) - int(quot(n))) % 2 == 0, range(0, N + 1))

M = N + 2                                                   # A379764 via the exact rewrite of the bilateral sum
Phi = [Fr(0)] * (M + 1); Phi[0] += 2; Phi[1] -= 1
for m in range(2, M + 2):
    base = [Fr(1)] + [Fr(0)] * M
    if m - 1 <= M: base[m - 1] -= 1
    P = power(base, m + 1, M)
    for i in range(M + 1 - (m + 1)): Phi[i + m + 1] += (-1) ** (m + 1) * P[i]
    e = m * (m - 1)
    if e <= M:
        b = [Fr(1)] + [Fr(0)] * M
        if m + 1 <= M: b[m + 1] -= 1
        P = power(inv(b, M), m - 1, M)
        for i in range(M + 1 - e): Phi[i + e] += P[i]
Q = [Fr(0)] * (N + 1)
for m in range(1, isqrt(N) + 1): Q[m * m] = (-1) ** (m + 1)
A = [Fr(0)] * (N + 1)
for _ in range(N + 1):
    new = Q[:]; Ak = A[:]
    for k in range(2, N + 1):
        Ak = mul(Ak, A, N)
        for i in range(N + 1): new[i] += 4 ** (k - 1) * Phi[k] * Ak[i]
    A = new
sq = lambda n: isqrt(n) ** 2 == n
report('a379764', A, 1, lambda n: (int(A[n]) % 2 == 1) == sq(n) and (((int(A[n]) - (-1) ** (n - 1)) % 4 == 0) == sq(n)) and (sq(n) or int(A[n]) % 4 == 0), range(1, N + 1))

def bil(F, p, n):                                           # A380681-A380686: the bilateral sum, exact rewrite
    F = (F + [Fr(0)] * (n + 1))[:n + 1]; tot = [Fr(0)] * (n + 1); tot[0] += 1
    for m in range(1, n + 1):
        sg = (-1) ** m
        if p * m <= n:
            xm = [Fr(0)] * (n + 1); xm[m] = Fr(1); t = power([F[i] + xm[i] for i in range(n + 1)], p * m, n)
            for i in range(n + 1 - p * m): tot[i + p * m] += sg * t[i]
        e = p * m * (m - 1)
        if e <= n:
            d = [Fr(1)] + [Fr(0)] * n
            for i in range(n + 1 - m): d[i + m] += F[i]
            t = power(inv(d, n), p * m, n)
            for i in range(n + 1 - e): tot[i + e] += sg * t[i]
    return tot
names = {2: 'a380681', 3: 'a380682', 5: 'a380683', 7: 'a380684', 11: 'a380685', 13: 'a380686'}
for p in (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15):
    n_ = 11 if p <= 7 else 8
    F = [Fr(1)] + [Fr(0)] * n_
    for k in range(1, n_ + 1):
        def val(G, k=k): t = bil(G, p, k + 1); return t[k + 1] - (p if k + 1 == 1 else 0)
        linsolve(F, k, val)
    full = bil(F, p, n_ + 1); full[1] -= p; eq = all(c == 0 for c in full)
    isprime = all(p % d for d in range(2, p))
    if isprime:
        ref = oeis(names[p], n_ + 1); m = [int(f) for f in F[:len(ref)]] == ref[:n_ + 1]
        print(f'p = {p:2d} ({names[p]}): equation {eq}, integral {integral(F)}, data match {m}'); ok &= eq and integral(F) and m
    else:
        q = min(d for d in range(2, p + 1) if p % d == 0); first = next(k for k in range(n_ + 1) if F[k].denominator != 1)
        pred = Fr(prod(1 + j * p for j in range(q)), factorial(q))
        good = eq and first == q - 1 and F[first] == pred
        print(f'p = {p:2d} composite: first non-integer coefficient index {first} = {F[first]} (predicted {pred}): {good}'); ok &= good
random.seed(7)
for _ in range(3):
    Y = [Fr(random.randint(-5, 5)) for _ in range(13)]
    z1 = all(c == 0 for c in bil(Y, 1, 12)); z3 = all(c.denominator == 1 and c.numerator % 3 == 0 for c in bil(Y, 3, 12))
    print(f'Psi_1(Y) = 0: {z1}; Psi_3(Y) = 0 mod 3: {z3}'); ok &= z1 and z3
print('RESULT', 'PASS' if ok else 'FAIL')
