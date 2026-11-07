# CRIM: four-part P-position candidate

`check.py` independently implements row deletion and conjugate-row-conjugate
column deletion from Bašić–Gottlieb–Krnc, *CRIM: A Natural Game on Integer
Partitions*, arXiv:2606.16828v1, Sections 2.1, 2.3 and 3.

Run with Python 3.10 or newer, using only the standard library:

```sh
python3 docs/reports/crim-question-one-four-parts/check.py
```

The output is stored in `check.out`. The exact integer evaluator reproduces
Proposition 10 on 465 two-part partitions. It then checks all 46,376 partitions
with at most four positive parts and width at most 30, including all 40,920
four-part partitions. It tests the candidate P-position criterion, equivalence
of source moves with column-height deletion, cell-count descent, and all thirteen
branches of an explicit winning-move strategy. The Grundy recursion does not
use the candidate criterion or the proposed strategy.

For a ≥ b ≥ c ≥ d ≥ 1 the candidate says P holds precisely when a is even and
at least one condition holds:

- d ≥ 2 and b,c are not both odd;
- d = c = 1 and either b is even or b = 1;
- d = 1, c ≥ 2 and b is odd;
- d = 1, c = 2 and either b ≥ 3 or a = 2.

This finite check supports a proposed answer to Section 5, Question 1. It is
neither an unbounded proof nor Lean kernel verification. The proposed settlement
requires induction on cell count over the full unbounded domain. The output also
records the boundary [1,1,1], which has Grundy value 1 and must be excluded from
the printed three-part characterization when it is used as an auxiliary result.

Source: https://arxiv.org/abs/2606.16828v1 . The program is an independent
implementation; no source code from the paper or thesis is included.
