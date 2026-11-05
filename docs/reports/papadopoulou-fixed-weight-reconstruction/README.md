# Fixed-weight reconstruction: majority succeeds while Algorithm 1 fails

`check.py` checks Conjecture 1 item 2 of Vivian Papadopoulou, Yonatan Yehezkeally,
V. Arvind Rameshwar and Antonia Wachter-Zeh, *When Majority-Vote Fails: Sequence
Reconstruction under Fixed-Weight Substitutions*, arXiv:2610.30001v1.
Source: <https://arxiv.org/abs/2610.30001v1>, Section IV-C.

Run with Python 3.9 or later, using only the standard library:

```sh
python3 docs/reports/papadopoulou-fixed-weight-reconstruction/check.py
```

The saved standard output is `output.txt`. Every assertion uses integer arithmetic
or finite set operations. The last line is `ALL_CHECKS_OK` on success.

The fidelity checks reproduce the source's seven-row, five-bit degeneracy example
and the complete two five-bit rows of its weighted-decoder table: 120 seven-view
sets with zero weighted failures and 90 ordinary-majority failures, and ten
nine-view sets with neither decoder failing.

The counterexample has 71 distinct nine-bit observations, all at distance four
from zero. Correct ordinary-majority column counts are
(28,29,33,33,32,33,32,32,32). The unique first weighted decision is coordinate
0 set to 0. Its budget update gives weight two to 28 rows and zero to the
other 43, making coordinate 1 the unique next maximizer with negative vote
−28; the other seven undecided votes are 20. The decoder sets that bit to 1.
No tie precedes this error. Full prescribed tie backtracking rejects every branch.
Enumeration of all 512 possible centers finds only zero.

The direct depth-first decoder counts 24 live states and 24 invalid child
choices. A separate implementation recomputes budgets from partial assignments,
merges identical states, and visits invalid terminals: its assignment DAG has
25 states, including eight invalid-budget terminals. These are different counting
conventions for the same rejecting computation.

Finite family checks cover zero padding to each length 9 through 14, enumerating
all possible centers at those lengths. All 512 XOR translations of the base
instance and four translations each at lengths 10 and 12 also fail while ordinary
majority is correct. These finite checks corroborate the padding and XOR argument;
they do not constitute an unbounded or Lean/kernel proof. The preregistration
route gives the general argument for every length at least nine and every
transmitted word. Item 1 of the conjecture is outside this refutation: the base
parameters satisfy N=71, N_wc=71, N_wmaj=101 and N_maj=113.
