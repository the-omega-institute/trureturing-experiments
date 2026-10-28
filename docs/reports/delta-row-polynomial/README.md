# DELTA row polynomial and Kurkov's accumulator loop

`python3 docs/reports/delta-row-polynomial/check.py` uses only Python's standard
library and integer arithmetic. Polynomials are integer coefficient tuples.
`output.txt` is the program's output.

The literal bivariate DELTA recurrence reproduces every displayed complete row
(n = 0 through 8) of [OEIS A111184](https://oeis.org/A111184/internal), revision
58. The program independently computes the sequential in-place algorithm in
Kurkov's Conjecture 3 for every pair 1 ≤ n ≤ m ≤ 20 and compares its result with
the DELTA polynomial. It also compares the quotient of rising-factorial series
H_(q+1)/H_q and return walks with transitions k → s ≤ k of weight 1 and
k → k+1 of weight q+k+1. All comparisons are identities in Z[q].

Additional checks evaluate q = -5,...,6, verify v_n(0) = n!, and recover the
repository's scalar recurrence U at q = 1. These finite results check source
fidelity, indexing, accumulator order, and the proposed proof's bridges. They
do not establish the unbounded conjecture.

Background: Paul Barry, arXiv:1804.06801v1, §4; Robert Cori,
arXiv:0812.0440v1, §1.2.2. The rising-factorial quotient is already recorded in
[OEIS A111528](https://oeis.org/A111528/internal). The target is the later
identification of the exact accumulator algorithm, rather than that quotient.
