# Barker's A336186 board-square conjectures

`check.py` uses Python 3 standard-library integer arithmetic, binomial coefficients,
and integer square roots. Run it with `python3 docs/reports/barker-a336186/check.py`.
Its complete output is in `check.out`; successful completion ends with `ALL_CHECKS_OK`.

The source is [OEIS A336186](https://oeis.org/A336186/internal), revision 38,
2025-08-20: the diagonal board, its example S(17) = 48841 = 221², and Colin Barker's
recurrence and ordinary generating-function conjectures. The program reproduces
the displayed board, both source boundary lists for the 17-by-17 block, the block
sum, and the ten published terms.

The proposed route reduces positive square-sum sizes to x² − 21y² = 15. It uses
the norm-preserving matrix T = ((55,252),(12,55)), its integral inverse, two seeds
(6,1) and (78,17), and positive descent for y ≥ 47. The checks verify the exact
matrix identities, all 46 base sizes, direct block sums for 0 ≤ m ≤ 64, and the
square/norm equivalence for every 1 ≤ m ≤ 2,000,000. They compare that entire
bounded enumeration with the two orbits, then test the recurrence at indices
5 through 200 and the formal generating-function identity at coefficients
0 through 200. Finite checks are evidence for the proposed route; they do not
prove the infinite classification, existence of the full source enumeration,
or either all-index conjecture. Those remain Lean proof obligations after
preregistration.
