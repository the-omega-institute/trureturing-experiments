"""Exact checks for Koutschan-Krattenthaler-Schlosser, arXiv:2401.08481v1, Conjecture 23, equation (10.12):
G_{2,1,2,0}(n) = 2 G_{7,4,7,-1}(n-1), with the binomial convention (1.2), together with the three lemmas of the
proposed proof (column generating functions, common multiplier h, triangular reduction)."""
from fractions import Fraction as F
import sympy as sp

def binom(x, p):
    if p < 0:
        return F(0)
    r = F(1)
    for s in range(p):
        r *= F(x - s)
    for s in range(1, p + 1):
        r /= s
    return r

def entry(a, b, g, d, i, j):
    return F(2) ** (i + b) * binom(i + 4 * j + g, 4 * j + a) + binom(-i + 4 * j + d, 4 * j + a)

def G(a, b, g, d, n):
    if n == 0:
        return F(1)
    M = sp.Matrix(n, n, lambda i, j: sp.Rational(entry(a, b, g, d, i, j)))
    return F(str(M.det()))

ok = True
for n in range(1, 13):
    lhs, rhs = G(2, 1, 2, 0, n), 2 * G(7, 4, 7, -1, n - 1)
    ok &= lhs == rhs
print("identity (10.12) for n = 1..12:", ok)

N = 48
t = sp.symbols('t')
def series_coeffs(expr, N):
    s = sp.series(expr, t, 0, N).removeO()
    return [F(str(sp.Poly(s, t).coeff_monomial(t ** k))) for k in range(N)]

lem1 = True
for j in range(5):
    aj = series_coeffs(2 / (1 - 2 * t) ** (4 * j + 3) + t ** (4 * j + 1) / (1 - t) ** (4 * j + 3), N)
    bj = series_coeffs(16 / (1 - 2 * t) ** (4 * j + 8) - t ** (4 * j) / (1 - t) ** (4 * j + 8), N)
    lem1 &= all(aj[i] == entry(2, 1, 2, 0, i, j) for i in range(N))
    lem1 &= all(bj[i] == entry(7, 4, 7, -1, i, j) for i in range(N))
print("Lemma 1 (column generating functions), j < 5, i <", N, ":", lem1)

X, Y, Z = sp.symbols('X Y Z')  # X = (1-2t)^(-4j), Y = (1-t)^(-4j), Z = t^(4j)
a = lambda k: 2 * X * (1 - 2 * t) ** (-4 * k - 3) + Z * t ** (4 * k + 1) * Y * (1 - t) ** (-4 * k - 3)
bj = 16 * X * (1 - 2 * t) ** (-8) - Z * Y * (1 - t) ** (-8)
h = 1 / ((1 - t) * (1 - 2 * t) * (1 - 2 * t + 2 * t ** 2))
lem2 = sp.simplify(h * (a(1) - a(0)) - t * bj) == 0
lem2 &= sp.expand(1 - (1 - 2 * t) ** 4 - 8 * t * (1 - t) * (1 - 2 * t + 2 * t ** 2)) == 0
lem2 &= sp.expand(t ** 4 - (1 - t) ** 4 + (1 - 2 * t) * (1 - 2 * t + 2 * t ** 2)) == 0
print("Lemma 2 (common multiplier, symbolic in j):", lem2)

hc = series_coeffs(h, 16)
lem3 = hc[0] == 1
for n in range(1, 16):
    A = sp.Matrix(n, n, lambda i, j: sp.Rational(entry(2, 1, 2, 0, i, j)))
    L = sp.Matrix(n, n, lambda i, k: sp.Rational(hc[i - k]) if k <= i else 0)
    R = sp.Matrix(n, n, lambda i, j: 1 if i == j else (-1 if j == i + 1 else 0))
    T = L * A * R
    lem3 &= T[0, 0] == 2 and all(T[0, j] == 0 for j in range(1, n))
    lem3 &= all(T[i, j] == sp.Rational(entry(7, 4, 7, -1, i - 1, j - 1)) for i in range(1, n) for j in range(1, n))
print("Lemma 3 (L A R has top row (2,0,...,0) and lower block B_{n-1}), n < 16:", lem3)

print("RESULT", "PASS" if (ok and lem1 and lem2 and lem3) else "FAIL")
