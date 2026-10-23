#!/usr/bin/env python3
"""Checks for D. Altizio, "Zeros of Stern polynomials in the complex plane", arXiv:2511.03847v1, Section 4.3, Remark 4.

(z)_n = 1 + z + ... + z^(n-1) for positive n; B = {|z - 2| <= 1}. Claim: |(z)_n| >= n|z| on B iff n in {5, 6, 7}.
On |z - 2| = 1 with u = |z|^2: z + conj(z) = (u + 3)/2 and |z - 1|^2 = (u - 1)/2, so P_n(u) = z^n + conj(z)^n satisfies
P_0 = 2, P_1 = (u + 3)/2, P_n = ((u + 3)/2) P_{n-1} - u P_{n-2}, and
D_n(u) = u^n + 1 - P_n(u) - n^2 u (u - 1)/2 = |z - 1|^2 (|(z)_n|^2 - n^2 |z|^2).
Exact (sympy): the source's factorizations (4.17) and (4.19) of Proposition 8; the certificates for n = 5, 6, 7 with
nonnegative coefficients in t = u - 1; the t^2 coefficient c_n = -n^2 (n - 5)(n - 7)/48 of D_n(1 + t) for n <= 60; and
D_n(1 + t_n) < 0 at t_n = 4^(-n-3) for every excluded n <= 60. Numeric: boundary and interior minima of
|(z)_n| - n|z| for n <= 20. Exits 1 on failure.
"""
import sys
import numpy as np
import sympy as sp

ok = True
u, t, r = sp.symbols('u t r')
N = 60
P = {0: sp.Integer(2), 1: (u + 3) / 2}
for n in range(2, N + 1):
    P[n] = sp.expand((u + 3) / 2 * P[n - 1] - u * P[n - 2])
D = lambda n: sp.expand(u ** n + 1 - P[n] - sp.Rational(n * n, 2) * u * (u - 1))
q = sp.symbols('q', real=True)
# the reduction identity at a symbolic boundary point z = 2 + e^{iq}
z = 2 + sp.exp(sp.I * q)
for n in (3, 5, 7):
    lhs = sp.simplify(sp.expand(D(n).subs(u, sp.expand(z * sp.conjugate(z)))) - sp.expand(sp.Abs(z - 1) ** 2 * (sp.Abs(sum(z ** k for k in range(n))) ** 2 - n * n * sp.Abs(z) ** 2)))
    good = abs(complex(lhs.subs(q, 0.37))) < 1e-9 and abs(complex(lhs.subs(q, 2.1))) < 1e-9
    print(f"reduction identity D_n(|z|^2) = |z - 1|^2 (|(z)_n|^2 - n^2 |z|^2) on |z - 2| = 1, n = {n}: {good}")
    ok &= good
f417 = sp.expand(sp.Rational(1, 32) * (31 * r**10 + 5 * r**8 + 10 * r**6 + 30 * r**4 + 135 * r**2 - 211) - sp.Rational(25, 2) * (r**4 - r**2) - (r**2 - 1)**3 * (31 * r**4 + 98 * r**2 + 211) / 32) == 0
f417b = sp.expand(u**5 + 1 - P[5] - sp.Rational(1, 32) * (31 * u**5 + 5 * u**4 + 10 * u**3 + 30 * u**2 + 135 * u - 211)) == 0
f419 = sp.expand(sp.Rational(1, 64) * (63 * r**12 + 6 * r**10 + 9 * r**8 + 20 * r**6 + 81 * r**4 + 486 * r**2 - 665) - 18 * (r**4 - r**2) - (r**2 - 1)**2 * (63 * r**8 + 132 * r**6 + 210 * r**4 + 308 * r**2 - 665) / 64) == 0
f419b = sp.expand(u**6 + 1 - P[6] - sp.Rational(1, 64) * (63 * u**6 + 6 * u**5 + 9 * u**4 + 20 * u**3 + 81 * u**2 + 486 * u - 665)) == 0
print(f"fidelity Proposition 8: (4.17) factorization {f417}, its left side equals |z^5 - 1|^2 {f417b}; (4.19) factorization {f419}, its left side equals |z^6 - 1|^2 {f419b}")
ok &= f417 and f417b and f419 and f419b
certs = {5: (32, t**3 * (340 + 160 * t + 31 * t**2)), 6: (64, t**2 * (48 + 1376 * t + 984 * t**2 + 384 * t**3 + 63 * t**4)), 7: (128, t**3 * (4704 + 4592 * t + 2716 * t**2 + 896 * t**3 + 127 * t**4))}
for n, (k, cert) in certs.items():
    g = sp.expand(k * D(n).subs(u, 1 + t) - cert) == 0
    print(f"certificate n = {n}: {k} D_{n}(1 + t) = {sp.factor(cert)}: {g}")
    ok &= g
cn_ok = all([sp.expand(D(n).subs(u, 1 + t)).coeff(t, i) for i in range(3)] == [0, 0, sp.Rational(-n * n * (n - 5) * (n - 7), 48)] for n in range(1, N + 1))
print(f"D_n(1) = D_n'(1) = 0 and c_n = -n^2 (n - 5)(n - 7)/48 for 1 <= n <= {N}: {cn_ok}")
ok &= cn_ok
excl = [n for n in range(1, N + 1) if n not in (5, 6, 7)]
w_ok = all(D(n).subs(u, 1 + sp.Rational(1, 4 ** (n + 3))) < 0 for n in excl)
print(f"D_n(1 + 4^(-n-3)) < 0 for every excluded n <= {N}: {w_ok}")
ok &= w_ok
th = np.linspace(0, 2 * np.pi, 200001)
rng = np.random.default_rng(20261023)
for n in range(1, 21):
    zb = 2 + np.exp(1j * th); mb = (np.abs(sum(zb ** k for k in range(n))) - n * np.abs(zb)).min()
    zi = 2 + np.sqrt(rng.random(200000)) * np.exp(2j * np.pi * rng.random(200000))
    mi = (np.abs(sum(zi ** k for k in range(n))) - n * np.abs(zi)).min()
    print(f"n = {n:2d}: min over |z - 2| = 1 of |(z)_n| - n|z| = {mb:+.3e}; over 2e5 interior samples {mi:+.3e}")
    ok &= (mb > -1e-9 and mi > 0) if n in (5, 6, 7) else (mb < 0)
print("ALL_CHECKS_OK" if ok else "CHECK_FAIL")
sys.exit(0 if ok else 1)
