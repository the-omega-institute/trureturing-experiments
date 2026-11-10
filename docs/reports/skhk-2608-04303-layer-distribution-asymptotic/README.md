# Random recursive tree layer distribution: exact checks

Run `python3 check.py` with Python 3.12 or later and SymPy 1.14.0. The checked output is in `check-output.txt`.

The program enumerates every labelled recursive-tree parent history at sizes 3 through 8 and reproduces both finite-size layer formulas of Theorem 7 in Stepanyants, Krapivsky, Hartle and Krioukov, arXiv:2608.04303v2. Leaf stripping is undirected and simultaneous; the terminal singleton is included, consistent with the source's layer-one fraction at size 3.

After substituting z=exp(-y), exact symbolic integrals reproduce B0=1, B1=1/2, B2=3-exp(1), and check the positive-product identity for individual masses at k=0,1. Exact rational Taylor coefficients check its differential residual for k=0 through 8 to degree 17. All evidence uses integers, fractions or exact symbolic expressions.

These checks do not prove the fixed-k limiting distribution, the identity for all k, the ABF branching-random-walk bounds, or Conjecture 1. The asymptotic formalization requires a research line, including those dependencies; a conditional weighted-integral theorem would not settle it.

Reference: H. Stepanyants, P. L. Krapivsky, H. Hartle and D. Krioukov, *The rank and layer distributions in random recursive trees*, arXiv:2608.04303v2, Section II and Section IV's finite-size proof. The positive-product identity is the proposed argument checked here, with recurrence/PWIT identification credited to Itoh (arXiv:1507.04379v6) and median/tail estimates to Addario-Berry and Ford (arXiv:1012.2544v3, Theorems 1.1–1.2).
