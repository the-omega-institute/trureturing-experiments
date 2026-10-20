"""Gawron-Ulas, arXiv:2202.05221v1, Conjecture 4.4: the second assertion fails for m = 3, k = 2.

S_{A,k}(n) = sum_i (-1)^i i^k c_A(i, n), with c_A(i, n) the number of A-partitions of n with i parts, and
U_3 = {1, 2, 6}. Checks, with exact integers:
1. the implementation reproduces the source's rational generating functions for k = 0, 1 (n < 80) and its sign
   statement for k = 1 (n >= 8);
2. S_{U_3,2}(12r + 2) = (-52 r^3 + 69 r^2 + 82 r + 9) / 3 for r <= 40, so (-1)^n S_{U_3,2}(n) < 0 at n = 12s + 38;
3. mechanism readings: for k <= 8 the sign of S_{U_3,k}(n) follows (-1)^floor((n+1)/2) on a window of large n.
"""
import sys
import sympy as sp


def S_direct(A, k, n):
    dp = [dict() for _ in range(n + 1)]
    dp[0][0] = 1
    for a in sorted(A):
        for m in range(a, n + 1):
            for i, c in list(dp[m - a].items()):
                dp[m][i + 1] = dp[m].get(i + 1, 0) + c
    return sum((-1) ** i * i ** k * c for i, c in dp[n].items())


def S126(n, k):
    t = 0
    for z in range(n // 6 + 1):
        for y in range((n - 6 * z) // 2 + 1):
            p = n - y - 5 * z
            t += (-1) ** p * p ** k
    return t


ok = True
A = [1, 2, 6]
x = sp.symbols('x')
F0 = 1 / ((1 + x) * (1 + x**2) * (1 + x**6))
F1 = -x * (3 * x**6 + 2 * x**5 - x**4 - x**3 + x**2 + x + 1) / ((1 + x)**2 * (1 + x**2) * (1 + x**6)**2)
N = 80
for k, F in ((0, F0), (1, F1)):
    ser = sp.Poly(sp.series(F, x, 0, N).removeO(), x)
    good = all(ser.coeff_monomial(x**n) == S_direct(A, k, n) for n in range(N))
    ok = ok and good
    print(f"k={k}: source generating function reproduced for n < {N}: {good}")
good = all((-1) ** ((n + 1) // 2) * S_direct(A, 1, n) > 0 for n in range(8, N))
ok = ok and good
print(f"k=1: source sign statement (-1)^floor((n+1)/2) S > 0 for 8 <= n < {N}: {good}")
for r in range(0, 41):
    n = 12 * r + 2
    val = S_direct(A, 2, n) if n <= 200 else S126(n, 2)
    form = sp.Rational(-52 * r**3 + 69 * r**2 + 82 * r + 9, 3)
    if val != form:
        ok = False
        print(f"closed form fails at r={r}")
    if r >= 3 and not (-1) ** n * val < 0:
        ok = False
    if r <= 6:
        print(f"n={n}: S_(U_3,2)(n) = {val}, closed form {form}, (-1)^n S < 0: {(-1) ** n * val < 0}")
print("closed form and violations checked for r <= 40")
for k, base in ((2, 400), (3, 400), (4, 800), (5, 1000), (6, 1500), (7, 3000), (8, 6000)):
    p4 = all((-1) ** ((n + 1) // 2) * S126(n, k) > 0 for n in range(base, base + 8))
    alt = sum(1 for n in range(base, base + 8) if (-1) ** n * S126(n, k) < 0)
    print(f"mechanism k={k}, n in [{base}, {base + 7}]: period-4 sign holds {p4}, (-1)^n violations {alt}/8")
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
