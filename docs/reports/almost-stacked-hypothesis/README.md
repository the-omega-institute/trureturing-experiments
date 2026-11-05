# Almost Stacked Hypothesis: exact-size threshold counterexample

Source: Tamás Csernák and Lajos Soukup, *The Almost Stacked Hypothesis:
A Conjectural Analogue of Sjöstrand's Cover Pebbling Theorem*,
arXiv:2610.09650v1, Hypothesis 3.2(1).
https://arxiv.org/html/2610.09650v1

Run `python3 check.py`. Python 3's standard library suffices; every calculation
uses integers. The deterministic output is stored in `check.out`.

On the connected simple graph with vertices 0 through 6 and edges
(0,1), (0,4), (0,5), (0,6), (1,2), (1,4), (2,4), (2,5), (3,4), (4,6), (5,6),
the least successful almost-stacked exact total is 12, whereas the least
successful unrestricted exact total is 15. In particular, all 448 almost-stacked
configurations of size 12 stack, but (0,0,0,9,0,3,0) does not. Its reachable set
has 225 states, is closed under every legal move, and contains no stacked state.
Almost-stacked failures return at totals 13 and 14.

The recursive search terminates because each move decreases total mass by one.
A separate bottom-up program enumerates all 170544 configurations of total at
most 15 and uses independently implemented adjacency-based transitions. It
agrees with the recursive evaluator on every almost-stacked configuration at
totals 2 through 15 and on the negative witnesses. Each of the 448 positive paths
is also checked move by move. Exhaustiveness is checked against all weak
compositions at each total. Theorem 4.10 for C4 at total 7 and Theorem 4.4 for P3
are reproduced as source-fidelity checks.

This refutes the printed least-exact-size assertion. It does not refute the
weaker assertion about almost-stacked obstructions at the final obstruction
level, the clearing clause, or the separate tree and cycle formulas. These are
computational results with a mathematical proof outline, not a Lean proof.
