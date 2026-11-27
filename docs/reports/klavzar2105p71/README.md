# Steiner general position monotonicity

Klavžar, Kuziak, Peterin and Yero ask about consecutive maximum Steiner general position numbers in Problem 7.1 of [arXiv:2105.08391v1](https://arxiv.org/abs/2105.08391v1).

Run `python3 check.py` with Python 3.12 or newer. The script uses only the standard library and exact integer arithmetic. The output from the recorded run is in `check.out`.

For each terminal set, the checker searches induced connected vertex sets in increasing cardinality. Every connected set of size `s` admits a spanning tree with `s - 1` edges. The checker takes the union over the entire minimum layer, so an extra selected vertex in any minimum tree rules out the selection. It then finds the largest selection avoiding the resulting forbidden sets by exhaustive search.

The fidelity cases apply Theorems 3.1 and 3.2 to paths, stars and cycles of orders 3 through 9. The script also checks the printed Corollary 4.4(v) at `r = 2, s = 5`, where the independent side gives the missing lower bound of 5 at `k = 3`. For the proposed family it checks both exact maxima, 2-connectivity, all 18 obstruction trees and a later increase for `2 ≤ m ≤ 8`. The uniform argument in the preregistration covers all `m ≥ 2`.
