"""Miska-Ulas, arXiv:2203.03942v1, Conjecture 4.8: a family X_n with sigma_3(X_n) = sigma_n(X_n) exceeding the bound Q(n).

For n >= 4, k = n - 4, Y = (1^k, 2, n - 1, c) with c = (3n^2 - 11n + 10)/2, d = sigma_3(Y) and X_n = (Y, d).
Checks, with exact integers, for 4 <= n <= 200: X_n is a nondecreasing positive tuple, sigma_3(X_n) = sigma_n(X_n),
x_n = D(n), m_n - Q(n) = (27/4) prod_{j=4}^{10} (n - j), and m_n > Q(n) exactly for n >= 11.
"""
import sys
from fractions import Fraction


def esym(xs):
    e = [1]
    for x in xs:
        e = [(e[i] if i < len(e) else 0) + (x * e[i - 1] if i >= 1 else 0) for i in range(len(e) + 1)]
    return e


def Q(n):
    return Fraction(3021 * n**6 - 77575 * n**5 + 920361 * n**4 - 6235705 * n**3 + 24764202 * n**2 - 53664304 * n
                    + 48986640, 12)


def D(n):
    return Fraction(27 * n**4 - 190 * n**3 + 471 * n**2 - 500 * n + 216, 12)


ok = True
for n in range(4, 201):
    k = n - 4
    assert (3 * n * n - 11 * n + 10) % 2 == 0
    c = (3 * n * n - 11 * n + 10) // 2
    Y = [1] * k + [2, n - 1, c]
    d = esym(Y)[3]
    X = Y + [d]
    E = esym(X)
    m = E[n]
    prod = Fraction(27, 4)
    for j in range(4, 11):
        prod *= (n - j)
    conds = [len(X) == n, min(X) >= 1, X == sorted(X), E[3] == E[n], d == D(n), m - Q(n) == prod,
             (m > Q(n)) == (n >= 11)]
    if not all(conds):
        ok = False
        print(f"n={n}: failed {conds}")
    if n <= 13 or n % 50 == 0:
        print(f"n={n}: sigma_3 = sigma_n = {m}, x_n = {d} = D(n), m - Q(n) = {m - Q(n)}")
print("X_11 =", [1] * 7 + [2, 10, 126, 16177])
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
