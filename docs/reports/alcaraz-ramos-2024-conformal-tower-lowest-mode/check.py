#!/usr/bin/env python3
"""Alcaraz-Ramos (arXiv:2402.09936v3; PRE 109 (2024) 044138), prediction (4.7) with (4.8), (4.10): for N = 3, lambda = 1,
M = 3m, M eps_i(M) -> pi 2 sqrt(3) (7/6 + i - 1), eps_i = z_i^(-1/3) over the positive roots of P_M = P_{M-1} - z P_{M-3}.

Scaled by u = m eps, the roots satisfy Q_m(u) = sum_k (-1)^k u^(3k) C(m+2+2k, 3k+2) / m^(3k+2) = 0, and Q_m -> F(u) =
sum_k (-1)^k u^(3k)/(3k+2)! = [e^(-u) + 2 e^(u/2) cos(sqrt3 u/2 - 2 pi/3)]/(3u^2). Checks (60-digit mpmath):
(1) the closed form of F equals the series at test points;
(2) for m = 50, 200, 800, 3200 the first zero of Q_m is a root of P_{3m} (sign change of the recurrence) and 3 u_1(m) increases
    towards 3 u_1(F) = 12.69962...; Richardson extrapolation from m = 800, 3200 agrees with it to 1e-3;
(3) 3 u_1(F) differs from the prediction 2 sqrt3 pi 7/6 = 12.69659... by more than 1e-3, and 3 u_1(F)/(2 sqrt3 pi) reproduces
    the source's Table I entry 1.1669439 to 1e-6;
(4) the first three zeros of F are not equally spaced.
Exit 0 and RESULT PASS when all hold.
"""
import mpmath as mp, sys
mp.mp.dps = 60
F = lambda u: (mp.e ** (-u) + 2 * mp.e ** (u / 2) * mp.cos(mp.sqrt(3) * u / 2 - 2 * mp.pi / 3)) / (3 * u ** 2)
G = lambda u: mp.nsum(lambda k: (-1) ** int(k) * u ** (3 * k) / mp.factorial(3 * k + 2), [0, mp.inf])
def Q(m, u): return mp.fsum((-1) ** k * u ** (3 * k) * mp.binomial(m + 2 + 2 * k, 3 * k + 2) / mp.mpf(m) ** (3 * k + 2) for k in range(m + 1))
def P(M, z):
    a = [mp.mpf(1)] * 3
    for n in range(1, M + 1): a.append(a[-1] - z * a[-3])
    return a[-1]
def bis(f, a, b, it=150):
    fa = f(a)
    for _ in range(it):
        c = (a + b) / 2; fc = f(c)
        if fa * fc <= 0: b = c
        else: a, fa = c, fc
    return (a + b) / 2
ok = True
e1 = max(abs(F(mp.mpf(u)) - G(mp.mpf(u))) for u in ('0.7', '3.1', '9.4')); ok &= e1 < mp.mpf('1e-30'); print('(1) closed form vs series max error', mp.nstr(e1, 3))
zs = []; u = mp.mpf('0.5'); f0 = F(u)
while len(zs) < 3:
    u2 = u + mp.mpf('0.01'); f1 = F(u2)
    if f0 * f1 < 0: zs.append(bis(F, u, u2))
    u, f0 = u2, f1
lim = 3 * zs[0]; pred = 2 * mp.sqrt(3) * mp.pi * mp.mpf(7) / 6
vals = []
for m in (50, 200, 800, 3200):
    u1 = bis(lambda x: Q(m, x), mp.mpf('3.9'), mp.mpf('4.6'))
    z = (m / u1) ** 3; sc = P(3 * m, z * (1 - mp.mpf('1e-20'))) * P(3 * m, z * (1 + mp.mpf('1e-20'))) < 0
    vals.append(3 * u1); print(f'(2) m={m}: 3 u_1(m) = {mp.nstr(3 * u1, 12)}, sign change of P_(3m) at z = (m/u1)^3: {sc}'); ok &= sc
rich = vals[-1] + (vals[-1] - vals[-2]) / 3; incr = all(vals[i] < vals[i + 1] for i in range(3))
print('(2) increasing:', incr, ' Richardson', mp.nstr(rich, 10), ' limit 3u_1(F)', mp.nstr(lim, 12)); ok &= incr and abs(rich - lim) < mp.mpf('1e-3')
gap = lim - pred; tab = lim / (2 * mp.sqrt(3) * mp.pi)
print('(3) limit - prediction =', mp.nstr(gap, 6), ' ratio', mp.nstr(tab, 10), '(Table I: 1.1669439)'); ok &= gap > mp.mpf('1e-3') and abs(tab - mp.mpf('1.1669439')) < mp.mpf('1e-6')
d1, d2 = zs[1] - zs[0], zs[2] - zs[1]; print('(4) spacings', mp.nstr(d1, 12), mp.nstr(d2, 12)); ok &= abs(d1 - d2) > mp.mpf('1e-4')
print('RESULT PASS' if ok else 'RESULT FAIL'); sys.exit(0 if ok else 1)
