# Exact checks for Gaunt's Bessel K tail conjecture

Source: Robert E. Gaunt, *Inequalities for some integrals involving modified
Bessel functions*, arXiv:1806.00524v2, Conjecture 2.4, Remark 2.5 and Appendix B;
Proc. Amer. Math. Soc. 147 (2019), 2937–2951, DOI: 10.1090/proc/14433.

Run `python3 check.py`. Python's standard library suffices; all arithmetic
uses integers and `fractions.Fraction`. The recorded output is in `check.out`.

The factorial coefficient formula for half-integer K functions independently
reproduces the source's three-term recurrence (20 cases) and weighted
derivative identity (42 cases, with positive and negative orders). For all
66 pairs a=A+1/2, b=B+1/2 with 2≤A≤12 and 0≤B≤A−2, set
ν=(a²−b²−1)/2 and n=a−ν. The source's parameter threshold is checked exactly.
The normalized strict derivative gap has zero constant and linear
coefficients and strictly positive remaining coefficients. Consequently
this finite family passes the derivative comparison for every x>0.

Additional checks cover rational parameter geometry (including b=0),
monomial secant/tangent comparisons, and the asymptotic coefficient matching
in Remark 2.5. The third normalized coefficient gap is
(ν²−b²)/2=(1−n)(2ν+n+1)/2>0.

These are finite algebra checks. The classical half-integer analytic
representation, convergence and differentiation of improper integrals,
infinite-series convexity and asymptotic remainder estimates are not
certified by this program. No general real-parameter proof or Lean proof is
claimed. The script uses no floating-point quadrature or approximations.
