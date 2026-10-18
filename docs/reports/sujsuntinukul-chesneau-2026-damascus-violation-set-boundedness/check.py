"""Sujsuntinukul–Chesneau arXiv:2601.00916v4, Section 5, question (ii): is I_n^m bounded for arbitrary m?
Classification: for m >= 2, I_n^m is unbounded iff (m-2) M_n >= 1, M_n = max_{x>0} (x^n-1)/(x^{n+1}+1);
explicitly bounded iff m <= 3, or m = 4 and n <= 4, or m in {5,6} and n = 1.
1. M_n via the root r_n of x^{n+1} - (n+1) x - n (M_n = n/((n+1) r_n)), compared with a dense grid maximum.
2. Lemma 5 exactly: M_1 = (sqrt2 - 1)/2 in (1/5, 1/4); M_2 = 1/3 at x = 2 with x^3 - 3x^2 + 4 = (x-2)^2 (x+1);
   the n = 4 identity x^5 - 2x^4 + 3 = (x - 8/5)^2 (x^3 + 6/5 x^2 + 32/25 x + 128/125) + 1183/3125; f_5(3/2) = 422/793 > 1/2.
3. Unbounded regimes: the family (t, r, ..., r, r^{-(m-2)}/t) has S > 0 at t = 10^3, 10^6 (exact Fractions for the
   rational choices r = 2 or 3/2 used in the issue; float for r = r_n), including the equality case (t,2,2,2,1/(8t)).
4. Bounded regimes: random product-one samples never give a violation with a coordinate above the bound B = K^(m-1).
5. Table of the classification for m <= 10, n <= 10.
"""
import math, random, sys
from fractions import Fraction as Fr
import sympy as sp

random.seed(20261017)
f = lambda n, x: (x ** n - 1) / (x ** (n + 1) + 1)

def r_root(n):
    lo, hi = 1.0, 3.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if mid ** (n + 1) - (n + 1) * mid - n < 0: lo = mid
        else: hi = mid
    return lo

def M(n): r = r_root(n); return n / ((n + 1) * r)

def bounded_formula(m, n): return m <= 3 or (m == 4 and n <= 4) or (m in (5, 6) and n == 1)

def main():
    fails = []
    for n in range(1, 21):
        grid = max(f(n, 1 + k / 20000) for k in range(1, 40000))
        if abs(grid - M(n)) > 1e-7: fails.append(("M", n))
    print("M_n for n = 1..8:", ", ".join(f"{M(n):.6f}" for n in range(1, 9)))
    x = sp.symbols('x')
    ok5 = (sp.simplify(sp.Rational(1, 2) * (sp.sqrt(2) - 1) - sp.Rational(1, 5)) > 0 and sp.Rational(1, 2) * (sp.sqrt(2) - 1) < sp.Rational(1, 4)
           and abs(M(1) - (math.sqrt(2) - 1) / 2) < 1e-12
           and sp.expand(x ** 3 - 3 * x ** 2 + 4 - (x - 2) ** 2 * (x + 1)) == 0 and abs(M(2) - 1 / 3) < 1e-12
           and sp.expand(x ** 5 - 2 * x ** 4 + 3 - ((x - sp.Rational(8, 5)) ** 2 * (x ** 3 + sp.Rational(6, 5) * x ** 2 + sp.Rational(32, 25) * x + sp.Rational(128, 125)) + sp.Rational(1183, 3125))) == 0
           and Fr(3, 2) ** 5 - 1 == Fr(211, 32) and (Fr(3, 2) ** 5 - 1) / (Fr(3, 2) ** 6 + 1) == Fr(422, 793) and Fr(422, 793) > Fr(1, 2))
    print(f"Lemma 5 exact identities: {ok5}")
    if not ok5: fails.append("lemma5")
    table = []
    for m in range(1, 11):
        row = []
        for n in range(1, 11):
            unb = m >= 2 and (m - 2) * M(n) >= 1 - 1e-15
            row.append('U' if unb else 'B')
            if unb == bounded_formula(m, n): fails.append(("table", m, n))
        table.append(f"m={m:2d}: " + ''.join(row))
    print("classification (n = 1..10, U unbounded, B bounded):"); print("\n".join(table))
    # exact equality-threshold family m=5, n=2
    for t in (Fr(10), Fr(1000), Fr(10 ** 6)):
        xs = [t, Fr(2), Fr(2), Fr(2), 1 / (8 * t)]
        assert math.prod(xs) == 1
        S = sum(f(2, v) for v in xs)
        if not S > 0: fails.append(("m5n2", t))
    # rational families: r = 2 for n = 1, m >= 7 and 2 <= n <= 4, m >= 5; r = 3/2 for n >= 5, m >= 4
    for m in range(4, 11):
        for n in range(1, 11):
            if bounded_formula(m, n): continue
            r = Fr(2) if n <= 4 else Fr(3, 2)
            for t in (Fr(10 ** 3), Fr(10 ** 6)):
                xs = [t] + [r] * (m - 2) + [1 / (r ** (m - 2) * t)]
                if math.prod(xs) != 1 or not sum(f(n, v) for v in xs) > 0: fails.append(("family", m, n, t))
    print("unbounded regimes (m <= 10, n <= 10): rational families give S > 0 at t = 10^3, 10^6:", not any(isinstance(z, tuple) and z[0] == "family" for z in fails))
    # bounded regimes: sampling
    worst = 0
    for m, n in ((3, 1), (3, 7), (4, 2), (4, 4), (5, 1), (6, 1)):
        delta = 1 - (m - 2) * M(n); K = 1 + 8 / delta; B = K ** (m - 1)
        for _ in range(60000):
            u = [random.gauss(0, random.choice((0.3, 1, 3, 8))) for _ in range(m - 1)]; u.append(-sum(u))
            xs = [math.exp(v) for v in u]
            if sum(f(n, v) for v in xs) > 0:
                worst = max(worst, max(xs) / B)
        if worst > 1: fails.append(("bounded", m, n))
    print(f"bounded regimes: largest violating coordinate / B over samples = {worst:.3e} (must be <= 1)")
    print("RESULT", "PASS" if not fails else f"FAIL {fails[:6]}")
    sys.exit(0 if not fails else 1)

if __name__ == "__main__":
    main()
