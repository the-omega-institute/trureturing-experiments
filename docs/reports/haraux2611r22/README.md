# Haraux sharp negative-potential threshold

Source: Alain Haraux, *An oscillation result for slightly perturbed trigonometric series*, arXiv:2611.04738v1, Remark 2.2. The proposed answer is the sharp bound c l² ≤ j², including equality, where j is the first positive zero of J₀.

Run `python3 check.py` with Python 3.12 or newer and SymPy 1.14.0. The stored output is `check.out`. The program uses integer, rational and symbolic exact arithmetic throughout.

The checks reproduce the source's c=3/4 obstruction; verify beta moments, transmutation coefficients and the weight-integral coefficients on the stated finite ranges; enclose the first Bessel zero using alternating series and a positive Bernstein polynomial; and certify an explicit two-mode witness at c=2/3. Rational Taylor bounds at 806 grid points and a derivative bound give trace >1/20 throughout [-π,π]. The PDE, boundaries and trace at x=π/√2 are checked symbolically. A symbolic cosh family checks the obstruction for c≥1.

These are exact finite checks and a finite-witness certificate. They do not prove the universal threshold, its endpoint, or the equivalence between the spectral and Sobolev energy formulations. The universal positive settlement is planned as an analytic argument; this experiment is supporting evidence.

The check program is original code under Apache-2.0. It copies no third-party program or Lean source.
