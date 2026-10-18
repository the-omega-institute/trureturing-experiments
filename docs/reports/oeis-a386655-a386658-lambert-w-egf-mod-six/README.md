# OEIS A386655–A386658: Hanna's mod 6 patterns for Σ (q^n x + LambertW(x))^n / n!

For q = 2, 3, 4, 5, P. D. Hanna conjectures that a(n) mod 6 repeats with period 6 from n = 6 on. The proof has two parts:
- **Closed form.** With e^(yW(x)) = Σ y(y−m)^(m−1) x^m/m! (Abel–Lagrange), a(n) = Σ_k C(n,k) q^(k(k+1)) (q^k − (n−k))^(n−k−1).
- **Periodicity.** Modulo a prime dividing q only the term k = 0 survives. Modulo a prime not dividing q the sum is a binomial transform of a period-6 function. Moreover (1+x)^9 = (1+x)^3 modulo (3, x^6 − 1), and (1+x)^2 = 0 modulo (2, x^2 − 1).

`check.py` (Python 3 with sympy) checks:
- the closed form against an exact series expansion for n ≤ N1;
- the published data;
- the four patterns for 6 ≤ n ≤ N2 (and already from n = 3);
- the two finite-ring identities.

Run `python3 check.py 12 400`. Exit 0 and `RESULT PASS` mean every check holds (about 9 minutes, mostly the series expansion).
