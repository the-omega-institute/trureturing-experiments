"""Checks for the factorial-boundary conjecture of Grau, Oller-Marcen and Varona, "A class of weighted Delannoy numbers",
arXiv:2501.09726v2 (Filomat 36(17) (2022)), Section 6: f(m,0) = m!, f(0,n) = n!,
f(m,n) = f(m-1,n) + f(m,n-1) + f(m-1,n-1); a_n = f(n,n); conjecture a_{n+1}/a_n - (n+1) -> 0.
Exact rational arithmetic; also the search seat's estimate a_n = 2 e^2 n! (1 + 2/n + O(n^-2))."""
import sys
from fractions import Fraction
from math import factorial, e

def table(N):
    f = [[0] * (N + 2) for _ in range(N + 2)]
    for m in range(N + 2): f[m][0] = factorial(m)
    for n in range(N + 2): f[0][n] = factorial(n)
    for m in range(1, N + 2):
        for n in range(1, N + 2):
            f[m][n] = f[m - 1][n] + f[m][n - 1] + f[m - 1][n - 1]
    return f

def main():
    N = 400
    f = table(N)
    print("first values a_n:", [f[n][n] for n in range(8)])
    print("source/OEIS check f(1,1), f(2,2), f(3,3):", f[1][1], f[2][2], f[3][3])
    prev = None
    for n in [5, 10, 20, 50, 100, 200, 400]:
        a, a1 = f[n][n], f[n + 1][n + 1]
        diff = Fraction(a1, a) - (n + 1)
        ratio = a / factorial(n) / (2 * e * e)
        corr = (ratio - 1) * n
        print(f"n={n}: a_(n+1)/a_n - (n+1) = {float(diff):.6e}; a_n/(2e^2 n!) = {ratio:.10f}; n*(that-1) = {corr:.6f}")
    diffs = [abs(float(Fraction(f[n + 1][n + 1], f[n][n]) - (n + 1))) for n in range(50, N)]
    mono = all(diffs[i + 1] <= diffs[i] for i in range(len(diffs) - 1))
    print("|difference| decreasing for 50 <= n < 400:", mono)
    print("EXIT=0")


def table_b(N, b):
    f = [[0] * (N + 2) for _ in range(N + 2)]
    for m in range(N + 2): f[m][0] = b(m)
    for n in range(N + 2): f[0][n] = b(n)
    for m in range(1, N + 2):
        for n in range(1, N + 2):
            f[m][n] = f[m - 1][n] + f[m][n - 1] + f[m - 1][n - 1]
    return f

def second():
    """The n^n-border conjecture of the same paragraph: f_{n+1,n+1}/f_{n,n} - e(n + 1/2) -> 0 (OEIS A346385)."""
    from decimal import Decimal, getcontext
    getcontext().prec = 60
    N = 300
    f = table_b(N, lambda n: n ** n)
    print("A346385 first values:", [f[n][n] for n in range(6)])
    E = Decimal(1).exp()
    for n in [5, 10, 20, 50, 100, 200, 300]:
        r = Decimal(f[n + 1][n + 1]) / Decimal(f[n][n])
        print(f"n={n}: a_(n+1)/a_n - e(n+1/2) = {float(r - E * (Decimal(n) + Decimal(1) / 2)):.6e}")

if __name__ == "__main__":
    main()
    second()
