# Oldekop Conjecture 1.8 exact check

`check.py` uses only integer matrices and SymPy exact determinants and ranks. It evaluates the proposed invariant

\[
\sum_{S\subseteq[d]}(-1)^{|S|+1}(\operatorname{rank} C[S]\bmod 2)
\]

for the all-one scaling and the paper's rank-three squared-difference scaling at `d = 3, …, 7`. The readings are `1, 4, 11, 26, 57` for the generic all-one matrices and `1, 3, 6, 10, 15` for rank three, matching the paper's generic value `2^(d−1)−d` and rank-three value `binom(d−1,2)`.

It then checks 32 exact diagonal-congruence variants of fixed zero-diagonal matrices in the rank classes used by the paper's Section 6 tables: the `d = 5` rank-four formula `5 + k`, and the `d = 6` rank-four, rank-five, and rank-six formulas `5 + k + r`, `4 + k + l + r`, and `5 + k + l + r`. Here `k` and `l` count nonzero principal minors of orders four and five, and `r` counts rank-three principal five-submatrices. The corresponding source table readings are d=4: 4 generically and 3 in the rank-three drop; d=5: 11 generically, 6 at rank three, and the rank-four/rank-five formulas; d=6: 26 generically and the rank-four/rank-five/rank-six formulas and lower bounds listed in Section 6. Every assertion is exact; no floating-point value is used as evidence.

Run with Python 3 and SymPy. The final output line is `ALL_CHECKS_OK`.
