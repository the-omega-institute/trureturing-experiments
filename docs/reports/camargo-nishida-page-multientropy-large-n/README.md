# Camargo–Nishida Page and multi-entropy large-N checks

Source: H. A. Camargo and M. Nishida, *Genuine Multi-Entropy of Fully Symmetric Gaussian States*, arXiv:2609.30754v2, §6 after Eqs. (63)–(64). The source requests analytic large-N formulas for the configurations in Eqs. (60)–(61).

Run `python3 check.py`. Python 3's standard library suffices. The stored output is `check.out`. All arithmetic uses integers and `fractions.Fraction`; no floating-point or numerical logarithm is evidence.

The program independently assembles the weighted replica-permutation matrices and evaluates their determinants by exact elimination. It reproduces the first three n=3 rows of Table 1 at four exact parameters each. It checks the proposed finite-N exponential entropy ratio on all ordered positive tripartitions with 3 ≤ N ≤ 9 at r = 1, 1/2, 1/3, 1/5 (336 cases), and cross-checks the reduction against the full 9N-dimensional Gaussian precision determinant for every ordered positive tripartition with N ≤ 5 at r = 1/2 (10 cases). Here r is the ratio of the constant-mode and transverse eigenvalues; the program verifies the source radical expression using its square and nonnegative sign. Thus every sample is a source state with a ≥ 1, with a determined exactly by N and r.

The claimed ratio is

\[
\exp(12\,\mathrm{GM}^{(3)}_3)=
\frac{[1+u(AB+AC+BC)]^2}
{[1+uA(B+C)][1+uB(A+C)][1+uC(A+B)]},\qquad
u=\frac{3(a^2-1)}{4(N-1)}.
\]

The checks also verify the Page/ME substitutions and algebraic factorizations yielding constants 400/27 and 27/2, and the bookkeeping in the proposed logarithmic error bounds.

These are finite exact checks and algebraic consistency checks. They do not prove the quantified large-N claim, the Gaussian integral bridge, or the logarithmic inequalities. The bridge is reused from the repository's existing Gaussian multi-entropy owner; no new Lean proof is part of this report.
