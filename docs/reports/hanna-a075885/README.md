# Hanna A075885 exact diagnostics

Run `python3 docs/reports/hanna-a075885/check.py` from the repository root.
The captured standard output is in `check.out`.

The fidelity check reproduces the 30 terms at indices 0 through 29 and the
example a(5) = 46 in [OEIS A075885](https://oeis.org/A075885/internal), revision
28 (2025-02-16 08:32:47), authored by Paul D. Hanna. The quoted numerical
sequence data are attributed to OEIS and its contributors; see the
[OEIS end-user license](https://oeis.org/LICENSE).

All arithmetic uses Python integers and `fractions.Fraction`. For every
1 ≤ n ≤ 1000 the program checks the terminal plateau, its geometric bound,
the sum/product bounds, exact floor multiplicities, and the floor-product
identities, including the rational-product identity recorded by Ridouane
Oudra in A010786. It also checks the finite summation-by-parts identity as
exact coefficients of formal logarithms and the rational floor-error and
telescoping-tail identities used in the proposed proof.

These are finite diagnostics for a proposed general proof. They establish
no infinite limit, literature novelty, or Lean proof. No floating point
numbers, numerical logarithms, or numerical roots are used.
