# Bangachev Open Problem 1: exact exponent check

Source: Kiril Bangachev, *A Short Note on Algorithms For The Noisy k-XOR
Problem*, arXiv:2610.27147v1, Problem 1, Lemma 1 and Open Problem 1.

Run `python3 docs/reports/bangache2610p1/check.py` from this repository.
The program uses only Python's standard library and exact integer and rational
arithmetic. Its complete output is in `check.out`; success ends with
`ALL_CHECKS_OK`.

The source fidelity check verifies the cycle-to-test signal cancellation in
4,840 cases and the exact Bernoulli XOR formula in 40 cases. The occupancy check
examines 6,684 ordered one-sparse matrices and 361,164 row subsets, verifies
`E[D] = choose(m,2)/n`, and checks that the maximum zero-sum subset size is
`2 sum floor(C_r/2) <= 2D`. Nine full H0/H1 enumerations independently check the
duplicate-pair tester's advantage formula.

The family uses exact source sizes `n=t^5`, `m=t^3` and positive noise
`eta=1/t^2`, with no rounding. The finite bounds give tester advantage at least
`1/24 > 1/100` for `t >= 2`, while Markov bounds the probability of any cycle
of length at least `n^(3/10)` by `1/sqrt(t)`. Rational exponent checks cover 90
parameter pairs in `1/2 < beta < 1`, `2 beta-1 < delta < beta`.

These finite checks support the stated asymptotic argument; they do not prove
its universal limit statements or certify any Lean theorem. The proposed
obstruction concerns the literal exact scale `n^(delta+o(1))`. It does not
refute an upper-bound reading: a duplicate pair is a cycle of length 2.
