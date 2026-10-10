# Zhang–Zhao: asymptotic normality of Speyer's g-polynomial coefficients for growing rank

`check.py` checks the finite ingredients of a proof route for Conjecture 4.8 of arXiv:2408.15506v1 (Acta Math. Sinica Engl. Ser. 42 (2026) 1087–1098): for every d = d(n) → ∞ with 1 ≤ d ≤ ⌊n/2⌋, the coefficients S_i(n,d) = (n−i−1)!/((d−i)!(n−d−i)!(i−1)!) of g_{n,d}(t) are asymptotically normal by local and central limit theorems.

- The identity g_{n,d}(t)/t = Σ_j C(d−1,j) C(n−d−1,j) (1+t)^j, coefficient by coefficient, for 4 ≤ n ≤ N and 2 ≤ d ≤ ⌊n/2⌋ (exact integers).
- The variance bound σ² ≥ (d−1)(n−d−1)/(4(n−2)) ≥ (d−1)/8 for the coefficient distribution (exact rationals), with the minimal ratio σ²/bound.
- For n ≤ 24: g_{n,d}(t)/t has exactly d−1 real zeros in (−∞, −1] counted by Sturm sequences, and does not vanish at −1 (sympy, exact).

Run: `python3 check.py 40` with sympy (exit 0; final line `pairs=361 N=40 bad=0 min_var_over_bound=1.026844`).
