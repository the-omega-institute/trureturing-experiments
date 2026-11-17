# Two-block depth-law counterexample check

This report checks Conjecture 15 of Noam Ben Shimon and Aryeh Lev
Zabokritskiy (Yohananov), *Extinction Depth and q-ary Error-Correcting Codes
for the Limited Permutation Channel*, arXiv:2607.19566v1, Section 3.3.

Run with Python 3.12 or later, without optimization:

```sh
python3 docs/reports/shimon2607c15/check.py
```

`check.py` uses exact strings and integers and has no external dependencies.
`check.out` contains its complete output, ending in `ALL_CHECKS_OK`.

The program reproduces the paper's first-extinction example and Proposition
12 for lengths 3 through 8, in both orientations. As a check of the channel
interpretation, it compares disjoint-swap balls with the source's literal
radius-one displacement-permutation balls for every binary word of lengths
1 through 6.

For the prefix-free pair 001 and 100000, it enumerates all legal block
prefixes and all disjoint adjacent-swap sets through depth 11. It finds
finite extinction depth 11, above the conjectured bound 8. It verifies the
explicit surviving pair at depth 10, the complete five-by-three prefix
sets at depth 11, and separation of every output by the count and position
of its ones. This is exact computational evidence, not a Lean proof.

The program is independently written from the mathematical definitions;
no third-party program or Lean source is copied or adapted. The paper's
companion program source was not read because an Apache-compatible code
license was not established.
