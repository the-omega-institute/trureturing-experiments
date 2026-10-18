"""OEIS A091190 / A091200 (comment): for x*A(x)^p = B(x*A(x^p)), B(x) = x/(1 - p*x), A(0) = 1,
"it appears that A(x) is an integer series only when p is prime".

Equivalently A(x)^p = A(x^p) / (1 - p x A(x^p)). The coefficient of x^n gives
p a_n = [x^n] (A(x^p) / (1 - p x A(x^p)) - A_<n(x)^p), defining A uniquely.
 - p prime: the bracket is 0 mod p (A_<n^p = A_<n(x^p) mod p), so A is integral.
 - p composite with least prime factor q: for n < p the right side is truncated 1/(1 - p x), so
   a_n = prod_{i<n} (1 + i p) / n!; the numerator of a_q is prime to q, so a_q is not an integer.
Checks: integrality to order 14 (p < 10) or 8 for every prime p <= 23; for every composite p <= 30 the
first non-integer coefficient is a_q and equals the closed form; the data of A091190 (p = 3) and A091200 (p = 5).
"""
from fractions import Fraction
from math import factorial, prod

def mul(a, b, n):
    c = [0] * (n + 1)
    for i, x in enumerate(a[:n + 1]):
        if x:
            for j, y in enumerate(b[:n + 1 - i]):
                if y: c[i + j] += x * y
    return c

def power(a, k, n):
    r = [1] + [0] * n
    for _ in range(k): r = mul(r, a, n)
    return r

def solve(p, N):
    a = [Fraction(1)] + [Fraction(0)] * N
    for n in range(1, N + 1):
        At = a[:n] + [Fraction(0)] * (N + 1 - n)
        Ap = power(At, p, n)
        sub = [Fraction(0)] * (n + 1)
        for k in range(0, n // p + 1): sub[k * p] = a[k]
        rhs = [Fraction(0)] * (n + 1); term = sub[:]
        for _ in range(n + 1):
            for i in range(n + 1): rhs[i] += term[i]
            term = mul(mul(term, [0, p] + [0] * (n - 1), n), sub, n)
        a[n] = (rhs[n] - Ap[n]) / p
    return a

def oeis(path, k):
    v = []
    for line in open(path):
        if line[:2] in ('%S', '%T', '%U'): v += [int(t) for t in line.split(None, 2)[2].strip().strip(',').split(',') if t]
    return v[:k]

ok = True
for p in range(2, 31):
    N = 14 if p < 10 else 8
    a = solve(p, N)
    isprime = all(p % d for d in range(2, p))
    integral = all(x.denominator == 1 for x in a)
    if isprime:
        print(f'p = {p:2d} prime: integral to order {N}: {integral}'); ok &= integral
    else:
        q = min(d for d in range(2, p + 1) if p % d == 0)
        first = next(n for n in range(N + 1) if a[n].denominator != 1)
        form = Fraction(prod(1 + i * p for i in range(q)), factorial(q))
        good = first == q and a[q] == form
        print(f'p = {p:2d} composite (q = {q}): first non-integer a_{first} = {a[first]}, closed form {form}: {good}'); ok &= good
for path, p in (('oeis-a091190.txt', 3), ('oeis-a091200.txt', 5)):
    ref = oeis(path, 14); a = solve(p, len(ref) - 1)
    m = [int(x) for x in a] == ref
    print(f'{path}: first {len(ref)} terms match: {m}'); ok &= m
print('RESULT', 'PASS' if ok else 'FAIL')
