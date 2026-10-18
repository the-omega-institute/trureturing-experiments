# Li–Morgan–Siegrist (J. Anal. 2026), Conjecture 6.1: the inflection point is not unique

Source: Y. Li, R. Morgan, M. Siegrist, "Concavity structure of fractional differential equations and zeros of Mittag–Leffler functions", The Journal of Analysis (2026), DOI 10.1007/s41478-026-01192-9 (open access, CC BY 4.0): Definition 1.1, Equations (6)–(7), Theorem 5.2, Conjecture 6.1.

Conjecture 6.1: for 1 < α < 2 and f as in Theorem 5.2, F″ has exactly one zero in (0, 1), and F″ changes sign there.

## Program

`check.py` (sympy and mpmath). Run `python3 check.py`; it prints `RESULT PASS` when every check below holds.

- **Exact curvature polynomial.** For α = 3/2 and a polynomial f, `Γ(3/2) x^{3/2} F″(x)` is an exact rational polynomial. The program computes the fractional integral of each monomial by the binomial expansion (no Gamma closed form) and checks it against the Beta formula.
- **Instance f = (1 − x)(1 + 5x¹²).** P(x) = 180640483/1764461790 + x/2 − x² + (10485760/676039)x¹³ − (4194304/260015)x¹⁴. It has signs +, −, +, − at 7/10, 3/4, 4/5 and 9/10, and three simple zeros in (0, 1), at 0.712340…, 0.780776… and 0.896775…. Mid-point quadrature of F″ at four points agrees to 1e−41.
- **Family f_A = (1 − x)(1 + Ax²⁰), A ∈ [19, 21].** P_A is affine in A and has signs +, −, +, − at 1/2, 3/4, 19/20 and 1 for A = 19 and A = 21.
- **Other orders.** For the same shape, the polynomial x^{3−α}F″ has three sign changes on a 1/4000 grid at α = 1.1 (k = 8, A = 10), α = 1.3 (k = 8, A = 5) and α = 1.7 (k = 30, A = 10). These cases are floating-point only.
- **Hypotheses of Theorem 5.2.** The program checks them for each forcing: f is a polynomial, f(1) = 0, f = (1 − x)q with q > 0 on [0, 1], and x³f′(x) → 0.

Reading: `run.log` (exit 0, `RESULT PASS`, 2.7 s).
