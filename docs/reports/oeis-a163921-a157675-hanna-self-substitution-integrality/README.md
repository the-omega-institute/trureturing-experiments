# OEIS A163921 and A157675: integrality of Hanna's self-substituting exponential

Check for the open-problem preregistration proving the conjecture of OEIS A163921 (P. D. Hanna) and the integrality conjecture of OEIS A157675.

Run: `python3 check.py [N]` (Python 3, standard library; exact rationals; default N = 16 coefficients).

For a sequence L, F(x) = exp(Σ_{n≥1} L(n) xⁿ/n) and G is the solution with G(0) = 1 of G(x) = exp(Σ_{n≥1} G(L(n)x)ⁿ xⁿ/n), computed by fixed-point iteration (coefficient n of the right side depends only on coefficients below n of G). The program checks:
- L = σ: F and G are integral through x^{N−1}, and G starts 1, 1, 2, 7, 41, 385, 5769 (A163921);
- L(n) = 2ⁿ: F = 1/(1 − 2x) and G are integral, and G starts 1, 1, 3, 19, 237, 5741, 270857, 25099497 (A157675);
- four random integral F = ∏(1 − x^d)^{−c_d} with c_d ∈ [−3, 3]: F and G are integral;
- the control L(n) = [n = 2] (F = exp(x²/2), not integral): G is not integral.

`run.log` holds the output of the run cited by the preregistration.
