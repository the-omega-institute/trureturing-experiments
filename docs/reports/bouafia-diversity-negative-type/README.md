# Negative-type diversity deletion

This report checks a proposed negative answer to question (2) in the unnumbered open problem immediately after Theorem 4 of Philippe Bouafia, *Magnitude and diversity of trees*, arXiv:2605.03681v1. The algorithm uses the intended nonnegative stopping condition from the theorem's proof. The printed reversed inequality in Algorithm 1 is not used.

Run `python3 check.py` with Python 3 and SymPy. All arithmetic is exact: Python integers, Fraction and SymPy rational functions. The recorded output is `output.txt`; a successful run exits 0 and ends with `ALL_CHECKS_OK`.

The fidelity check reproduces Proposition 5's weighting and magnitude formulas for a four-leaf weighted star, and verifies Theorem 4's returned minimizer with the optimality criterion of Proposition 3. The main check builds the five binary vectors, verifies their Hamming metric and squared-Euclidean embedding identities, runs the deletion rule with an independent Fraction Gaussian solver, and checks the exact competing probability at q = 19/20 against SymPy. It also verifies the parametric solve, determinant and second-order residual and negative Schur-complement entry identities and exact Bernstein positivity certificates on 19/20 ≤ q ≤ 1. The algorithm and strict improvement conclusions use q < 1; q = 1 is used only to bound polynomials and take limits. All 31 nonempty principal minors at q = 19/20 are positive.

For each 19/20 ≤ q < 1, put c = -log(q) and d = cH. The binary-vector squared-distance identity gives the square-root embedding sqrt(c)x, and exp(-cH) = q^H. These real analytic implications are explained mathematically by those identities; Python does not numerically approximate exp, log or sqrt, and the program is not a Lean proof. The main Lean result and its negation have not been proved in this report.

The interval failure is specific to the positive-scale exponential kernel. The same configuration has an optimal scale-zero distance-matrix output. Irreversible deletion and re-entry as a general failure mechanism already appear in Huntsman, arXiv:2509.19565v3, Appendix G. This report does not claim that mechanism as new. The Euclidean-distance question (1), separately preregistered in main-repository issue #25735, is outside this report's settlement claim.

Source: https://arxiv.org/abs/2605.03681v1 (CC BY 4.0). The check program and output are an independent calculation; no source-paper code was copied.
