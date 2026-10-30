# Mandal–Manna distance signless Laplacian / Dilworth checks

Source: S. Mandal and P. Manna, *The Distance Laplacian and Distance Signless
Laplacian Spectra of C-Graphs*, [arXiv:2610.08599v1](https://arxiv.org/abs/2610.08599v1),
Example 5.4 and Conjecture 5.5. The graph is reconstructed from the paper's
complement/disjoint-union recurrence; distances are computed independently by BFS.

Run `python3 check.py`. Python 3's standard library suffices. All numerical
operations use integers or `fractions.Fraction`; no floating point is used.
The expected output is stored in `output.txt`, ending in `ALL_CHECKS_OK`.

The fidelity check reproduces Example 5.4's quotient matrix, its exact
characteristic polynomial, the full matrix's algebraic multiplicity two at 8,
the excluded set {5,7}, and Dilworth number two (two chains plus an incomparable
pair). The decimal roots in the paper are replaced by their exact quadratic.

For each integer r from 1 through 12, the candidate sequence is (r,1,...,1),
with 2r+5 trailing ones and order 3r+5. The program checks the independent
isolated/dominating construction, connectedness, every distance, transmission
formulas, two integer eigenvectors at 4r+4, a nonzero independence minor,
exact rational nullity, every ordered pair of the proposed vicinal chain,
and exclusion from the source's local-value set. It also checks that merging
the initial true-twin class gives a different local eigenvalue. At r=1 it
reproduces the matrix, vectors, transmissions, edge count, and factored
characteristic polynomial of the eight-vertex counterexample.

These are finite exact checks, not a Lean proof or an exhaustive literature
screen. The infinite family requires the separately stated symbolic argument.
