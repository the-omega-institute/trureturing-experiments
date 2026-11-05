# Singleton-loop nSSP classification checks

Source: Sara Koljančić and Polona Oblak, *Combinatorial aspects of the non-symmetric Strong Spectral Property for graphs*, arXiv:2604.21552v1, Problem 1.6 and the singleton subproblem following Example 5.10; ELA 42, 755–779, DOI 10.13001/ela.2026.10523.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The expected final line is `ALL_CHECKS_OK`; the complete exact readings are in [check-output.txt](check-output.txt).

The fidelity checks reproduce the displayed forbidden pair and the allowing adjacency matrix in Example 5.9, and the full-column-rank verification for the displayed nonsymmetric matrix in Example 5.10. The verifier uses an arbitrary real test matrix X, with transpose in the commutator and no extra diagonal condition.

For every predicted allowed pair (n,k), 1 ≤ n ≤ 12, the program chooses positive integer edge products with coprime branch continuants. It computes the exact determinant of the n-parameter characteristic-coefficient Jacobian and checks the characteristic polynomial directly. All 63 determinants are nonzero. For orders at most seven, a separate arbitrary-X linear rank calculation also checks nSSP directly. The rational nonsymmetric witness uses upper edges 1 and lower edges equal to the selected positive products.

The odd-order/even-loop obstruction is checked on all 15 such pairs through order eleven with nonzero signed directed edge weights and explicit left and right null vectors. The failure-of-requirement construction is checked for all 25 pairs of branch lengths r,s in {2,...,6}, using positive symmetric branches sharing eigenvalue 1, exact radicals, and the alternating-sign commutant. These matrices have the exact connected singleton-loop support.

Only exact integers, rational numbers, polynomials, and algebraic radicals are used. These finite checks test the proposed general argument; they do not prove the all-orders trichotomy or universal requirement at the boundary positions. The report contains independent check code and result data, not a Lean proof. The check program is original repository experiment code; source attribution is given above.
