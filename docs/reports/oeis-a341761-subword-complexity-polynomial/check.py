"""OEIS A341761 (S. Guo): S(n, x) = sum over words w of length n over an x-letter alphabet of the number of distinct
nonempty factors (contiguous subwords) of w; the rows are the coefficients of the polynomial S(n, x).
Conjecture: the coefficient of x^n in S(n, x) is n(n+1)/2.

Proof checked here: grouping words by the set of letters used, S(n, x) = sum_{j=1..n} C(x, j) * F(n, j), where F(n, j)
is the total factor count over words using exactly j given letters; C(x, j) has degree j, so only j = n contributes to
x^n, with C(x, n) = x^n/n! + ... and F(n, n) = n! * n(n+1)/2 (all n letters distinct: every factor is distinct).
The script computes S(n, x) by brute force for x <= n + 1, interpolates the polynomial, compares with the published
rows, and checks the leading coefficient and F(n, n), for n <= N.
"""
import sys
from itertools import product
from fractions import Fraction
from math import comb, factorial
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
def factors(w):
    return len({w[i:j] for i in range(len(w)) for j in range(i + 1, len(w) + 1)})
def S(n, x):
    return sum(factors(w) for w in product(range(x), repeat=n))
def interp(points):  # Lagrange -> coefficient list
    n = len(points); coeffs = [Fraction(0)] * n
    for i, (xi, yi) in enumerate(points):
        basis = [Fraction(1)]; denom = Fraction(1)
        for j, (xj, _) in enumerate(points):
            if j == i: continue
            basis = [Fraction(0)] + basis
            for k in range(len(basis) - 1): basis[k] -= xj * basis[k + 1]
            denom *= (xi - xj)
        for k in range(n): coeffs[k] += yi * basis[k] / denom
    return coeffs
DATA = [[0, 1], [0, -1, 3], [0, 0, -3, 6], [0, -1, 1, -6, 10], [0, 2, -6, 4, -10, 15]]
ok = True; rows = []
for n in range(1, N + 1):
    pts = [(x, S(n, x)) for x in range(0, n + 2)]
    c = interp(pts)
    row = [int(v) for v in c[:n + 1]]
    rows.append(row)
    if any(v.denominator != 1 for v in c) or c[n + 1:] != [0] * (len(c) - n - 1): ok = False
    if row[n] != n * (n + 1) // 2: ok = False
    distinct = sum(factors(w) for w in product(range(n), repeat=n) if len(set(w)) == n)
    if distinct != factorial(n) * n * (n + 1) // 2: ok = False
ok_data = rows[:len(DATA)] == DATA
print(f"rows n<={N}: {rows}")
print(f"matches data: {ok_data}; degree n with leading coefficient n(n+1)/2 and F(n,n) = n! n(n+1)/2: {ok}")
good = ok and ok_data
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
