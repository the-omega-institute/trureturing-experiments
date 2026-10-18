"""OEIS A380066 (Hanna): 1/2 = Sum_{n in Z} x^n A^n (A^n + x)^(n-1) (A + x^n)^(n-1), a(0) = 1. Conjecture: a(n) = 2 mod 4, n > 0.
This is the case c = 1 of the frozen D5.S1.Recurrence.Bilateral.BilateralProductModFour family (exponent c n - 1),
whose theorems assume c >= 2. Route (GPT Pro r386q, checked): the k = -m summand is (xA)^(m^2) /
((1 + x A^m)^(m+1) (1 + x^m A)^(m+1)), so the equation is 1/2 = 1/((1+x)(1+A)) + S(A) with S(A) integral and
contracting; at the fixed point (1+x)(1+A) = 2/(1 - 2S) = 2 (mod 4), so A = (1 - x)/(1 + x) (mod 4).
Checks: the rewrite reproduces A380066 and (c = 2) A381362; integrality and a(n) = 2 mod 4 for n <= 18;
A = (1 - x)/(1 + x) mod 4 coefficientwise.
"""
from fractions import Fraction as Fr

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

def residual(A, c, n):
    x = [Fr(0), Fr(1)] + [Fr(0)] * (n - 1); xA = mul(x, A, n)
    tot = mul(inv([Fr(1), Fr(1)] + [Fr(0)] * (n - 1), n), inv([1 + A[0]] + A[1:], n), n)
    for m in range(1, n + 1):
        Am = power(A, m, n); xm = [Fr(0)] * (n + 1); xm[m] = Fr(1)
        t1 = [a + (1 if i == 1 else 0) for i, a in enumerate(Am)]
        t2 = [xm[i] + A[i] for i in range(n + 1)]
        P = mul(power(xA, m, n), mul(power(t1, c * m - 1, n), power(t2, c * m - 1, n), n), n)
        d1 = [Fr(1)] + [Fr(0)] * n
        for i, v in enumerate(mul(x, Am, n)): d1[i] += v
        d2 = [Fr(1)] + [Fr(0)] * n
        for i, v in enumerate(mul(xm, A, n)): d2[i] += v
        Nm = mul(power(xA, c * m * m, n), inv(mul(power(d1, c * m + 1, n), power(d2, c * m + 1, n), n), n), n)
        tot = [a + b + d for a, b, d in zip(tot, P, Nm)]
    tot[0] -= Fr(1, 2)
    return tot

def solve(c, N):
    A = [Fr(1)] + [Fr(0)] * N
    for n in range(1, N + 1):
        A[n] = Fr(0); e0 = residual(A, c, n)[n]; A[n] = Fr(1); e1 = residual(A, c, n)[n]
        A[n] = -e0 / (e1 - e0)
    return A

def oeis(path, k):
    v = []
    for line in open(path):
        if line[:2] in ('%S', '%T', '%U'): v += [int(t) for t in line.split(None, 2)[2].strip().strip(',').split(',') if t]
    return v[:k]

ok = True
N = 18
A = solve(1, N)
ref = oeis('oeis-a380066.txt', 10)
m = [int(a) for a in A[:10]] == ref; print('A380066: first 10 terms match:', m); ok &= m
integral = all(a.denominator == 1 for a in A); print('integral to order', N, ':', integral); ok &= integral
mod4 = all(int(a) % 4 == 2 for a in A[1:]); print('a(n) = 2 mod 4 for 1 <= n <=', N, ':', mod4); ok &= mod4
geo = all((int(A[n]) - (2 * (-1) ** n if n else 1)) % 4 == 0 for n in range(N + 1))
print('A = (1 - x)/(1 + x) (mod 4) coefficientwise:', geo); ok &= geo
B = solve(2, 10)
m2 = [int(b) for b in B[:10]] == [1, 2, 46, 862, 20414, 526106, 14519710, 419293310, 12527971550, 384222183226]
print('c = 2 reproduces A381362 (frozen family):', m2); ok &= m2
print('RESULT', 'PASS' if ok else 'FAIL')
