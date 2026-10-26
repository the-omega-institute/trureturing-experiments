"""Gaunt (Proc. AMS 147 (2019) 2937-2951, arXiv:1806.00524v2) Conjecture 2.4: numerical check.

For n < 1, nu >= 1 - n + sqrt(2(1-n)), alpha = nu + n - sqrt((nu+n)^2 - 2nu - 1), q = nu + n - alpha:
  (2.17)  int_x^inf K_{nu+n}(t) t^{-nu} dt > K_q(x) x^{-nu}   for all x > 0;
  optimality: for r = |nu + n - beta| > q the inequality with K_r fails for large x (tested r = q + 1, x up to 400); for r < q, K_r < K_q;
  equality at n = 1.
Usage: PYTHONPATH=<mpmath> python3 gaunt-c24.py   (exit 0 iff all readings hold)
"""
import sys
from mpmath import mp, besselk, quad, inf, sqrt, mpf
mp.dps = 50
def T(nu, p, x): return quad(lambda t: besselk(p, t) * t**(-nu), [x, x + 1, x + 10, inf])
ok = True
grid = [mpf('0.01'), mpf('0.1'), mpf('0.5'), 1, 2, 5, 10, 20]
for n in [mpf('0.5'), 0, mpf('-1'), mpf('-3')]:
    nmin = 1 - n + sqrt(2 * (1 - n))
    for nu in [nmin, nmin + mpf('0.3'), nmin + 2]:
        p = nu + n; q = sqrt(max(0, p**2 - 2 * nu - 1))
        gaps = [T(nu, p, x) - besselk(q, x) * x**(-nu) for x in grid]
        good = all(g > 0 for g in gaps)
        ok &= good
        # optimality: r slightly above q fails at large x
        r = q + 1; xs = [50, 100, 200, 400]
        rel = [(T(nu, p, x) - besselk(r, x) * x**(-nu)) / (besselk(r, x) * x**(-nu)) for x in xs]
        ok &= rel[-1] < 0
        lower = besselk(q + mpf('-0.5') if q >= mpf('0.5') else 0, 1) <= besselk(q, 1)
        ok &= lower
        print(f"n={float(n):+.2f} nu={float(nu):.4f} q={float(q):.4f} min_rel_gap={float(min(g/(besselk(q,x)*x**(-nu)) for g,x in zip(gaps,grid))):.3e} r=q+1 rel_gap at x=50..400: {[f'{float(v):+.2e}' for v in rel]}")
# equality at n = 1
for nu in [mpf('-0.7'), 0, mpf('0.5'), 3]:
    e = max(abs(T(nu, nu + 1, x) - besselk(nu, x) * x**(-nu)) for x in [mpf('0.3'), 1, 4])
    print(f"n=1 nu={float(nu)} max|T - x^-nu K_nu| = {float(e):.2e}")
    ok &= e < mpf('1e-15')
sys.exit(0 if ok else 1)
