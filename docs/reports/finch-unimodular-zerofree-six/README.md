# Finch unimodular zerofree order-six counterexample

Source: Steven Finch, *Small Matrices with Small Inverses: Unimodular Zerofree Cases*, arXiv:2605.03691v4, 12 June 2026, Sections 1 and 5. The source defines both alpha and beta as maximum absolute entries and requires the matrix and its inverse to have no zero entries. Conjecture 4 excludes order-six matrices with alpha = beta = 3.

Run `python3 docs/reports/finch-unimodular-zerofree-six/check.py` with Python 3. The script uses standard-library integer arithmetic and `fractions.Fraction`, with no floating point or external dependencies. Its stored stdout is `output.txt`; successful completion ends with `ALL_CHECKS_OK`.

The fidelity check exhausts all nonzero-entry order-two integer matrices with maximum absolute entries bounded by 2 or 3. Exact determinants and inverses, followed by independent signed row and column permutations with the source's structural ordering, reproduce one orbit at alpha = beta = 2 and three at alpha = beta = 3, including exactly the source's representatives.

For the order-six candidate, the program checks AB = BA = I, nonzero entries, both maxima equal to 3, and determinants -1. Determinants are computed independently by fraction-free Bareiss elimination and all 720 Leibniz terms. Rational Gauss-Jordan elimination independently recovers each displayed inverse.

The family is A = D + k uu^T and B = D - k ss^T, where s has equally many +1 and -1 entries, D = diag(s), and u consists of ones. The program checks 72 pairs with m = 1 through 8 and k = 2 through 10, as well as all balanced sign placements for m = 1 through 4 and k = 2 through 4 (294 further pairs). For each m = 1 through 8, it also checks every polynomial coefficient in both inverse identities over the integer polynomial ring. Empty dimension, k = 1 and an unbalanced signature supply boundary controls.

These are exact finite computations. They certify the displayed finite counterexample by the implemented algorithms, but are not a Lean kernel proof or a proof for every dimension and parameter. The symbolic family argument belongs in the preregistration and subsequent formalization. The program does not claim a new general rank-one inverse formula or exhaustive literature priority.
