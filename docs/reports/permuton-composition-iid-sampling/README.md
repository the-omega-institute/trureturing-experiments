# Permuton composition under independent iid sampling

Exact finite checks for the iid-sampling question in Section 11 of Christian
Gorski, *On the Composition of Permutons*, arXiv:2610.12945v1.

Run `python3 check.py`. The program uses only integers and `fractions.Fraction`;
there are no external dependencies. The recorded output is `check.out`.

The fidelity checks reproduce the source's compatibility of checkerboard
composition with permutation composition at a common size, the uniform law of
a product of independent uniform permutations, and absorption by the uniform
permuton. All 15,017 permutation pairs of sizes 1 through 5 are checked.

The sampling check verifies that matching the second-coordinate ranks of Q to
the first-coordinate ranks of P gives the induced permutation sigma composed
with tau. It checks all coordinate-order types and both sample relabelings for
sizes 1 through 3 (1,313 cases).

Exact order-statistic expectations are computed in two independent ways: by
binomial tails and by integration of the beta-density polynomial. The diagonal
and antidiagonal examples include singular permutons. For the diagonal pair,
the expected matched lower-half rectangle mass also agrees with a separate
calculation using two independent binomial threshold counts; the finite checks
satisfy squared error at most 1/(8n). A noncommuting pair of three-cell rational
stationary kernels distinguishes the required Q-then-P rectangle mass 1/18
from the reverse mass 1/6. The sample sizes checked are 1, 2, 4, 8, 16, 32, 64.

These computations test source fidelity, composition order, the rank-matching
identity, and finite consequences of the proposed settlement. They establish
no almost-sure or asymptotic convergence theorem. A general measure-theoretic
argument is required for the proposed positive settlement.
