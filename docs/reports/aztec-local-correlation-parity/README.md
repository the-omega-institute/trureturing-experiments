# Aztec local-correlation parity check

This exact-arithmetic check independently reconstructs the Kasteleyn determinant certificate for the Aztec diamond orders `n = 4, 5, 6, 7` in arXiv:2610.22974v1. It computes the full determinant and the four-cell-deleted minors for the seven radius-3 orbit representatives, forms the uniform horizontal-pair probabilities, and checks the published central value `P_H(4) = 9/128`, the four published contrasts, and the predicted alternating signs.

Run from this directory with `python3 check.py`. SymPy evaluates the determinants over Gaussian integers; all displayed comparisons use integers and `fractions.Fraction`, with no floating-point evidence.
