# Mourad–Abbas Conjecture 5.1: exact checks

Source: Bassam Mourad and Hassan Abbas, *On the symmetric doubly stochastic matrices that are determined by their spectra*, arXiv:1310.1273v1, Sections 1, 2 and 5. The conjecture concerns uniqueness up to simultaneous row/column permutation among similar symmetric doubly stochastic matrices of positive trace.

Run `python3 docs/reports/mourad-abbas-c51/check.py` with Python 3 and SymPy 1.14.0. All evidence uses integers, rational numbers, or exact symbolic algebra. The stored stdout is `output.txt`; success exits 0 and ends with `ALL_CHECKS_OK`.

The fidelity checks reproduce the stated spectrum of C_n for orders 2 through 8 and the source's cospectral non-permutation-similar direct-sum example for k = 3 through 6. The candidate checks cover 25 rational instances of A = direct_sum_r(t I_k + (1-t) C_k), with k = 4 through 7 and r = 2 or 3, including k = 4, r = 2, t = 1/40. They verify stochasticity, trace, the two-eigenvalue polynomial, the shifted idempotent, all admissible component-size partitions, and symmetric four-cycle perturbations. Symbolic scalar identities check the general formulas and the perturbation on the ray from the identity. Boundary controls at t = 1/k² exhibit unequal-block cospectral mates for k = 4 through 6.

These checks support a preregistration draft. They do not prove spectral determination against every real matrix, replace the general proof outline, or establish a Lean settlement. The planned argument uses the classical nonnegative-idempotent classification (Peter Flor, *On groups of non-negative matrices*, Compositio Mathematica 21 (1969), 376–382, Theorem 2), followed by the component-size constraint and segment exclusions. No source text or third-party implementation is included.
