"""OEIS A207969/A207970 and A208034 (Hanna 2012): integrality of exp(sum c_n x^n/n) for c_n = 5 F_n^(2k), 2 P_n^(2k).

Integrality of exp(sum c_n x^n/n) is equivalent to the local Dold congruences
c_(p m) = c_m (mod p^(v_p(m)+1)) for all primes p and m >= 1.
Route: 5^k F_n^(2k) = (phi^n - phi'^n)^(2k) and 8^k P_n^(2k) = (alpha^n - beta^n)^(2k) are integer combinations of
traces of powers of integer matrices (Dold); dividing by 5^(k-1) (resp. 2^(3k-1)) is harmless away from 5 (resp. 2),
and at 5 (resp. 2) both sides of each congruence are divisible by the required power since v_5(F_n) = v_5(n)
(resp. v_2(P_n) = v_2(n)). Odd powers and the missing factor fail at an explicit prime.
"""
from fractions import Fraction

def expser(c, N):
    a = [Fraction(1)] + [Fraction(0)] * N
    for n in range(1, N + 1):
        a[n] = sum(Fraction(c[k]) * a[n - k] for k in range(1, n + 1)) / n
    return a

def seq(r, N):
    s = [0, 1]
    while len(s) <= N: s.append(r * s[-1] + s[-2])
    return s

def v(p, n):
    e = 0
    while n % p == 0: n //= p; e += 1
    return e

def dold_ok(c, primes, M):
    for p in primes:
        for m in range(1, M + 1):
            if p * m >= len(c): break
            if (c[p * m] - c[m]) % p ** (v(p, m) + 1): return (p, m)
    return None

N = 80; F = seq(1, 2 * N); P = seq(2, 2 * N)
integral = lambda a: all(x.denominator == 1 for x in a)
ok = True
print('Fibonacci, conjecture (5 F_n^(2k)):')
for k in range(0, 9):
    c = [5 * F[n] ** (2 * k) for n in range(2 * N + 1)]
    i = integral(expser(c, N)); d = dold_ok(c, [2, 3, 5, 7, 11, 13], 40)
    print(f'  k = {k}: integral to order {N}: {i}; local Dold congruences (p <= 13, m <= 40): {d is None}'); ok &= i and d is None
print('Fibonacci, notes (expected to fail):')
for k in range(0, 5):
    c = [5 * F[n] ** (2 * k + 1) for n in range(2 * N + 1)]
    w = (c[3] - c[1]) % 3; print(f'  5 F^(2k+1), k = {k}: integral {integral(expser(c, 40))}; c_3 - c_1 mod 3 = {w}'); ok &= w != 0 and not integral(expser(c, 40))
for k in range(1, 5):
    c = [F[n] ** (2 * k) for n in range(2 * N + 1)]
    w = (c[5] - c[1]) % 5; print(f'  F^(2k), k = {k}: integral {integral(expser(c, 40))}; c_5 - c_1 mod 5 = {w}'); ok &= w != 0
print('Pell, conjectures (1)-(3), k >= 1:')
for k in range(1, 9):
    c1 = [2 * P[n] ** (2 * k) for n in range(2 * N + 1)]
    c2 = [2 * P[n] ** (2 * k - 1) for n in range(2 * N + 1)]
    c3 = [P[n] ** (2 * k) for n in range(2 * N + 1)]
    i1 = integral(expser(c1, N)); d1 = dold_ok(c1, [2, 3, 5, 7, 11, 13], 40)
    w2 = (c2[3] - c2[1]) % 3; w3 = (c3[2] - c3[1]) % 2
    i2 = integral(expser(c2, 40)); i3 = integral(expser(c3, 40))
    print(f'  k = {k}: (1) integral {i1}, Dold {d1 is None}; (2) integral {i2}, c_3 - c_1 mod 3 = {w2}; (3) integral {i3}, c_2 - c_1 mod 2 = {w3}')
    ok &= i1 and d1 is None and not i2 and w2 != 0 and not i3 and w3 != 0
print('valuations: v_5(F_n) = v_5(n) and v_2(P_n) = v_2(n) for n <= 160:',
      all(v(5, F[n]) == v(5, n) and v(2, P[n]) == v(2, n) for n in range(1, 161)))
ok &= all(v(5, F[n]) == v(5, n) and v(2, P[n]) == v(2, n) for n in range(1, 161))
def oeis(path, k):
    vals = []
    for line in open(path):
        if line[:2] in ('%S', '%T', '%U'): vals += [int(t) for t in line.split(None, 2)[2].strip().strip(',').split(',') if t]
    return vals[:k]
for path, c in (('oeis-a207969.txt', [5 * F[n] ** 4 for n in range(30)]), ('oeis-a207970.txt', [5 * F[n] ** 6 for n in range(30)]),
                ('oeis-a208034.txt', [2 * P[n] ** 2 for n in range(30)])):
    ref = oeis(path, 12); m = [int(x) for x in expser(c, 11)] == ref
    print(f'{path}: first 12 terms match: {m}'); ok &= m
print('RESULT', 'PASS' if ok else 'FAIL')
