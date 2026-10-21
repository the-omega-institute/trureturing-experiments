"""Kagey, arXiv:2112.05281v1, Section 5.2: Conjectures 5.1 (inversions) and 5.2 (major index), exact checks.

Brute force over S_n for n <= 9:
- Conjecture 5.1: E[pi(1) | inv = k] = M(n+1, k)/M(n, k) for every k < n <= 8 (both sides exact rationals);
- Conjecture 5.2: the count A_k(n) and first-letter sum B_k(n) of {maj = k} agree, for every k <= 8 and k <= n <= 9 with a
  nonempty conditioning set, with the search seat's inclusion-exclusion polynomials, and Q_k = k! A_k and
  P_k = (k+1)! B_k are monic with integer coefficients and degrees k and k+1; the source's displayed cases (5.3)-(5.6)
  are reproduced.
"""
import itertools
import sys
from fractions import Fraction
from math import factorial, prod
import sympy as sp

x = sp.symbols("x")
N = 9
maj_stats, inv_stats = {}, {}
for n in range(1, N + 1):
    for p in itertools.permutations(range(1, n + 1)):
        maj = sum(i + 1 for i in range(n - 1) if p[i] > p[i + 1])
        inv = sum(1 for i in range(n) for j in range(i + 1, n) if p[i] > p[j])
        for table, key in ((maj_stats, (n, maj)), (inv_stats, (n, inv))):
            s = table.setdefault(key, [0, 0])
            s[0] += 1
            s[1] += p[0]

ok = True
for k in range(0, 8):
    cases = [n for n in range(k + 1, N) if (n, k) in inv_stats and (n + 1, k) in inv_stats]
    good = all(Fraction(inv_stats[(n, k)][1], inv_stats[(n, k)][0]) == Fraction(inv_stats[(n + 1, k)][0], inv_stats[(n, k)][0])
               for n in cases)
    ok = ok and good
    print(f"Conjecture 5.1, k={k}: E = M(n+1,k)/M(n,k) for n in {cases}: {good}")


def falling(m):
    return sp.prod([x - i for i in range(m)])


def C_U(U):
    if not U:
        return sp.Integer(1)
    U = sorted(U)
    d = [U[0]] + [U[j] - U[j - 1] for j in range(1, len(U))]
    return falling(U[-1]) / prod(factorial(t) for t in d)


def D_U(U):
    return sp.Integer(1) if not U else (x + 1) * C_U(U) / (min(U) + 1)


displayed = {1: (x**2 + x - 2, x - 1), 2: (x**3 - x - 6, x**2 - x - 2), 3: (x**4 + 6 * x**3 - 13 * x**2 - 18 * x, x**3 - 7 * x),
             4: (x**5 + 20 * x**4 - 45 * x**3 - 80 * x**2 - 16 * x, x**4 + 2 * x**3 - 13 * x**2 - 14 * x)}
for k in range(1, 9):
    A = B = sp.Integer(0)
    for r in range(1, k + 1):
        for S in itertools.combinations(range(1, k + 1), r):
            if sum(S) != k:
                continue
            for t in range(len(S) + 1):
                for U in itertools.combinations(S, t):
                    sign = (-1) ** (len(S) - len(U))
                    A += sign * C_U(U)
                    B += sign * D_U(U)
    Q, P = sp.Poly(sp.expand(factorial(k) * A), x), sp.Poly(sp.expand(factorial(k + 1) * B), x)
    shape = Q.LC() == 1 and P.LC() == 1 and Q.degree() == k and P.degree() == k + 1 and \
        all(c.is_integer for c in Q.all_coeffs() + P.all_coeffs())
    data = [n for n in range(k, N + 1) if (n, k) in maj_stats]
    match = all(A.subs(x, n) == maj_stats[(n, k)][0] and B.subs(x, n) == maj_stats[(n, k)][1] for n in data)
    shown = k not in displayed or (sp.expand(P.as_expr() - displayed[k][0]) == 0 and sp.expand(Q.as_expr() - displayed[k][1]) == 0)
    ok = ok and shape and match and shown
    print(f"Conjecture 5.2, k={k}: monic integer Q (deg {Q.degree()}) and P (deg {P.degree()}): {shape}; "
          f"matches brute force at n in {data}: {match}; source display reproduced: {shown}")
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
