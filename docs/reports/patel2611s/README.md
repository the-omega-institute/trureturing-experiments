# Patel movement-budget dynamic regret checks

The program uses only Python standard-library integers and exact rational
arithmetic. Run `python3 docs/reports/patel2611s/check.py`; the expected output
is stored in `check.out`, ending in `ALL_CHECKS_OK`.

Source: Ridham Patel, “When Is Tracking Affordable? A Movement-Budget
Feasibility Threshold for Dynamic Regret”, arXiv:2611.02812v1, Appendix F
(the OGD estimate) and Appendix B (the open question following Theorem 2).

The fidelity check reproduces the source's exact OGD substitutions, including
the step-size threshold and the movement bound. The generic finite sign
calculation cancels the coefficients of arbitrary history-dependent
block-start coordinates for 1 through 12 signs. The constructed-policy checks
enumerate every sign pattern for the displayed small phase/coordinate/block
parameters and verify the block lower bound, zero-loss common minimizers,
Euclidean-domain feasibility and hard movement feasibility. Actions change
one coordinate at a time, so their movement is an exact absolute difference.
Coordinates are normalized by a = D/(2 sqrt(m)); radicals cancel in these
checks. The parameter check covers T=1 through 128, a quarter-grid of P/D,
and eight budgets S/D, verifying the floor and ceiling choices and factor
256 using certified rational bounds on square roots.

These checks support the proposed general argument. They do not enumerate all
learner strategies, prove the minimax theorem, or constitute a Lean proof.
The intended theorem concerns the fixed Euclidean ball in dimension T and
comparators that minimize every current loss. It leaves the minimax rate in
each fixed small dimension unresolved.
