"""Xu-Yang (J. Inequal. Appl. 2026:98) Conjecture 7.1: numerical sign check.

R(x) = e^{x^2/2} int_x^inf e^{-t^2/2} dt (Mills ratio), Y_{a,b} = (x^2+a)R^2 + bxR,
u_n(a,b) from (3.7), D_n(x,a) = Y_{a,a-2}(x) - (a - 1 + sum_{k<=n} (-1)^k u_k / x^{2k}).
Conjecture: (-1)^n D_n(x,-1) is CM on (0,inf), i.e. (-1)^{n+m} D_n^{(m)} > 0.
Expected reading: holds for n = 0 and all n >= 2 (n <= 6 tested, m <= 4), fails for n = 1 (already m = 0).
Usage: PYTHONPATH=<mpmath> python3 xu-yang-c71.py  (exit 0 iff the expected readings hold)
"""
import sys
from mpmath import mp, mpf, erfc, sqrt, exp, pi, gamma, diff, factorial
mp.dps = 60
def R(x): return sqrt(pi/2) * exp(x**2/2) * erfc(x/sqrt(2))
def u(n, a, b):
    s = sum(2**k * gamma(k + mpf(1)/2) / (gamma(mpf(1)/2) * factorial(k)) for k in range(n))
    return (b + 1) * 2**n * gamma(n + mpf(1)/2) / gamma(mpf(1)/2) + (n - a) * factorial(n - 1) * s
a = -1; b = a - 2
assert u(1, a, b) == b + 2 - a and u(2, a, b) == 3*b + 7 - 2*a
print("u_1..u_6(-1,-3):", [int(u(k, a, b)) for k in range(1, 7)])
def Dn(n, x):
    return (x**2 + a) * R(x)**2 + b * x * R(x) - (a - 1 + sum((-1)**k * u(k, a, b) / x**(2*k) for k in range(1, n + 1)))
xs = [mpf('0.3'), mpf('0.7'), 1, 2, 4, 8, 15]
ok = True
for n in range(0, 7):
    worst = None; signs_ok = True
    for x in xs:
        for m in range(0, 5):
            v = (-1)**(n + m) * diff(lambda t: Dn(n, t), x, m)
            if v <= 0: signs_ok = False
            worst = v if worst is None or v < worst else worst
    expect = (n != 1)
    print(f"n={n}: all (-1)^(n+m) D_n^(m) > 0 on grid (m<=4): {signs_ok}  (expected {expect}); min value {float(worst):.3e}")
    ok &= signs_ok == expect
sys.exit(0 if ok else 1)
