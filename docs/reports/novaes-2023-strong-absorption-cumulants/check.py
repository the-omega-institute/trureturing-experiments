#!/usr/bin/env python3
"""Checks for Novaes, arXiv:2305.18235v1 (J. Stat. Phys. 190, 182 (2023)), Section 5, conjecture (58):
(-1)^n k_n(gamma) = (n-1)!/(M^(2n-2) gamma^(2n)) - (2n-1) n (n-1)!/(M^(2n-2) gamma^(2n+1)) + O(gamma^(-2n-2)).

Model (Grabsch, arXiv:1909.01002v2, Eqs. (3)-(4), beta = 2, N = M): with Gamma_i = gamma + x_i the joint law of (x, y)
has density proportional to exp(-M(sum x + sum y)) Delta(x)^2 Delta(y)^2 prod_{i,j} (gamma + x_i - y_j) on x_i >= 0,
0 <= y_j <= gamma, and tau_W = M^-1 sum_i 1/(gamma + x_i) in dwell-time units.

M = 1, high-precision quadrature:
1. the moment E[tau_W^2] agrees with Novaes Eq. (53) up to O(gamma^-6);
2. for n = 2..5 the cumulant error against (58), times gamma^(2n+2), stays bounded as gamma doubles (50..800);
3. for n = 3 that scaled error approaches the gamma^-8 coefficient -6(41 M^2 + 10)/M^6 = -306 of Novaes Eq. (56).
Prints ALL_OK on success.
"""
from mpmath import mp, quad, exp, inf, mpf, factorial
mp.dps = 50

def moments(g, kmax):
    w = lambda x, y: exp(-x - y) * (g + x - y)
    inner = lambda x: quad(lambda y: w(x, y), [0, g])
    Z = quad(inner, [0, inf])
    return [quad(lambda x: inner(x) / (g + x) ** k, [0, inf]) / Z for k in range(1, kmax + 1)]

def cumulants(m):  # m[0] = E X, ..., from raw moments by the standard recursion
    n = len(m); mm = [mpf(1)] + list(m); k = [mpf(0)] * (n + 1)
    from math import comb
    for r in range(1, n + 1):
        k[r] = mm[r] - sum(comb(r - 1, j - 1) * k[j] * mm[r - j] for j in range(1, r))
    return k

ok = True
gammas = [50, 100, 200, 400, 800]
scaled = {n: [] for n in range(2, 6)}
for g0 in gammas:
    g = mpf(g0); m = moments(g, 5); k = cumulants(m)
    p2 = 1/g**2 - 2/g**3 + 4/g**4 - 12/g**5
    e53 = float((m[1] - p2) * g**6)
    for n in range(2, 6):
        pred = factorial(n - 1) / g**(2*n) - (2*n - 1) * n * factorial(n - 1) / g**(2*n + 1)
        scaled[n].append(float(((-1)**n * k[n] - pred) * g**(2*n + 2)))
    print('gamma=%d  Eq.(53) err*g^6=%.4f  ' % (g0, e53) + '  '.join('n=%d: %.3f' % (n, scaled[n][-1]) for n in range(2, 6)))
    if abs(e53) > 100: ok = False
for n in range(2, 6):
    s = scaled[n]
    if not (abs(s[-1]) < 2 * abs(s[-2]) + 1 and abs(s[-1] - s[-2]) < abs(s[-2] - s[-3]) + 1e-9):
        print('NOT_CONVERGING n=%d' % n, s); ok = False
# n = 3: (-1)^3 k_3 = ... so k_3's gamma^-8 coefficient -306 appears with sign +306 in (-1)^3 k_3 - pred
if abs(scaled[3][-1] - 306) > 6:
    print('EQ56_MISMATCH', scaled[3][-1]); ok = False
print('ALL_OK' if ok else 'FAIL')
